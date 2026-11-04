"""
Query SigCom LINCS (L1000 chemical perturbations) with an up/down gene signature
and aggregate the per-signature hits into one score per drug.

Usage
-----
  python sigcom_query.py                                # runs the built-in dummy muscle signature
  python sigcom_query.py --up up.txt --down down.txt --name muscle_EE_4h
  python sigcom_query.py --up up.txt --down down.txt --cell-lines SKL ASC   # tissue-matched check

up.txt / down.txt: one gene symbol per line.
Outputs (in --outdir): <name>_signatures.csv (one row per LINCS signature hit)
                       <name>_drugs.csv      (one row per drug, aggregated across cell lines)

Endpoints and field names follow the SigCom LINCS tutorial and Consensus Appyter;
check the printed sanity output on the first run (see SIGN CONVENTION below).
"""
import argparse
import os
import sys

import pandas as pd
import requests

METADATA_API = "https://maayanlab.cloud/sigcom-lincs/metadata-api"
DATA_API = "https://maayanlab.cloud/sigcom-lincs/data-api/api/v1"
DATABASE = "l1000_cp"  # chemical perturbations

# Illustrative acute-endurance-exercise muscle response (~1-4 h post). NOT real MoTrPAC output:
# use it only to test the pipeline. It is CREB/cAMP-heavy, so cAMP agonists (forskolin/colforsin,
# PDE inhibitors, beta-adrenergic agonists) should come out as mimickers. Glucocorticoid targets
# (PER1, DUSP1, NFIL3) are included too, so dexamethasone-like drugs may also show up.
DUMMY_UP = """NR4A3 NR4A1 NR4A2 PPARGC1A ATF3 FOS FOSB JUNB EGR1 HES1 MYC MAFF SIK1 CREM
CCN1 CCN2 VEGFA HMOX1 ANKRD1 PDK4 HK2 PFKFB3 DUSP1 ZFP36 IER3 GADD45G CDKN1A HSPA1A
HSPA1B DNAJB1 SERPINE1 SOCS3 MT2A TNFRSF12A ADM NFIL3 PER1 KLF4""".split()
DUMMY_DOWN = """TXNIP ARRDC2 ARRDC3 ARRDC4 MSTN SREBF1 INSIG1 SCD FASN ELOVL6""".split()


def read_genes(path):
    with open(path) as f:
        return [g.strip() for g in f if g.strip()]


def resolve_genes(symbols):
    """Map gene symbols to SigCom entity UUIDs. Returns {symbol: uuid}."""
    payload = {"filter": {"where": {"meta.symbol": {"inq": symbols}}, "fields": ["id", "meta.symbol"]}}
    res = requests.post(f"{METADATA_API}/entities/find", json=payload, timeout=60)
    res.raise_for_status()
    return {e["meta"]["symbol"]: e["id"] for e in res.json()}


def rank_two_sided(up_ids, down_ids, limit):
    payload = {"up_entities": up_ids, "down_entities": down_ids, "limit": limit, "database": DATABASE}
    res = requests.post(f"{DATA_API}/enrich/ranktwosided", json=payload, timeout=300)
    res.raise_for_status()
    return res.json()["results"]


def signature_metadata(uuids, chunk=500):
    rows = []
    fields = ["id", "meta.pert_name", "meta.cell_line", "meta.pert_dose", "meta.pert_time", "meta.local_id"]
    for i in range(0, len(uuids), chunk):
        payload = {"filter": {"where": {"id": {"inq": uuids[i:i + chunk]}}, "fields": fields}}
        res = requests.post(f"{METADATA_API}/signatures/find", json=payload, timeout=120)
        res.raise_for_status()
        for s in res.json():
            m = s.get("meta", {})
            rows.append({"uuid": s["id"], "drug": m.get("pert_name"), "cell_line": m.get("cell_line"),
                         "dose": m.get("pert_dose"), "time": m.get("pert_time"), "local_id": m.get("local_id")})
    return pd.DataFrame(rows)


def tidy_hits(results, meta):
    df = pd.DataFrame(results).merge(meta, on="uuid", how="left")
    # SIGN CONVENTION: following the SigCom Consensus Appyter, z_sum = z-up + z-down,
    # and z_sum > 0 means a mimicker (the drug moves genes the same way as the query).
    # Check this on the first run: with the dummy signature, forskolin/colforsin
    # should land on the mimicker side. If they show up as reversers, flip SIGN.
    SIGN = 1
    df["z_sum"] = SIGN * (df["z-up"] + df["z-down"])
    df["direction"] = df["z_sum"].apply(lambda z: "mimicker" if z > 0 else "reverser")
    return df


def aggregate_by_drug(df, min_cell_lines=2):
    """One score per drug, robust to cell-line specificity.

    Per drug x cell line, take the strongest signal (across doses/times).
    Across cell lines: median score, number of cell lines hit, and sign consistency.
    Drugs absent from the top-`limit` results count as 'no hit', not zero, so
    n_cell_lines also reflects how widely a drug was tested. Cross-check the
    ranking against the Enrichr consensus library, which has one entry per drug.
    """
    per_cl = (df.loc[df.groupby(["drug", "cell_line"])["z_sum"].apply(lambda s: s.abs().idxmax())]
                [["drug", "cell_line", "z_sum"]])
    agg = per_cl.groupby("drug").agg(
        n_cell_lines=("cell_line", "nunique"),
        median_z=("z_sum", "median"),
        frac_mimic=("z_sum", lambda s: (s > 0).mean()),
        cell_lines=("cell_line", lambda s: ",".join(sorted(map(str, s)))),
    ).reset_index()
    agg["consistency"] = agg["frac_mimic"].apply(lambda f: max(f, 1 - f))
    agg["call"] = agg["median_z"].apply(lambda z: "mimicker" if z > 0 else "reverser")
    agg = agg[agg["n_cell_lines"] >= min_cell_lines]
    return agg.sort_values("median_z", ascending=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--up"); ap.add_argument("--down")
    ap.add_argument("--name", default="dummy_muscle_EE")
    ap.add_argument("--limit", type=int, default=1000, help="number of top signatures to return")
    ap.add_argument("--cell-lines", nargs="*", help="keep only these cell lines (e.g. SKL ASC)")
    ap.add_argument("--min-cell-lines", type=int, default=2)
    ap.add_argument("--controls", nargs="*", default=["forskolin", "colforsin", "rolipram", "isoproterenol",
                                                      "clenbuterol", "metformin", "GW-501516", "simvastatin",
                                                      "atorvastatin", "dexamethasone"])
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()

    up = read_genes(a.up) if a.up else DUMMY_UP
    down = read_genes(a.down) if a.down else DUMMY_DOWN
    mapping = resolve_genes(up + down)
    missing = [g for g in up + down if g not in mapping]
    print(f"Mapped {len(mapping)}/{len(up) + len(down)} genes. Unmapped: {missing or 'none'}")
    up_ids = [mapping[g] for g in up if g in mapping]
    down_ids = [mapping[g] for g in down if g in mapping]
    if len(up_ids) < 5 or len(down_ids) < 5:
        sys.exit("Need at least 5 mapped genes per direction.")

    results = rank_two_sided(up_ids, down_ids, a.limit)
    meta = signature_metadata([r["uuid"] for r in results])
    hits = tidy_hits(results, meta)
    if a.cell_lines:
        hits = hits[hits["cell_line"].isin(a.cell_lines)]
        print(f"Restricted to cell lines {a.cell_lines}: {len(hits)} signatures")

    drugs = aggregate_by_drug(hits, a.min_cell_lines if not a.cell_lines else 1)
    os.makedirs(a.outdir, exist_ok=True)
    hits.to_csv(os.path.join(a.outdir, f"{a.name}_signatures.csv"), index=False)
    drugs.to_csv(os.path.join(a.outdir, f"{a.name}_drugs.csv"), index=False)

    print(f"\n{len(hits)} signatures, {hits['drug'].nunique()} drugs, {hits['cell_line'].nunique()} cell lines")
    print("Cell lines seen:", ", ".join(sorted(map(str, hits["cell_line"].dropna().unique()))))
    cols = ["drug", "n_cell_lines", "median_z", "consistency"]
    print("\nTop mimickers:\n", drugs.head(15)[cols].to_string(index=False))
    print("\nTop reversers:\n", drugs.tail(15)[::-1][cols].to_string(index=False))
    ctrl = drugs[drugs["drug"].str.lower().isin([c.lower() for c in a.controls])]
    print("\nControls found:\n", ctrl[cols + ["call"]].to_string(index=False) if len(ctrl) else "none in results")


if __name__ == "__main__":
    main()