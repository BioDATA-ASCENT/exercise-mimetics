#!/usr/bin/env python3
"""Export FDR-significant human muscle-cell signatures from selected GEO studies."""

from __future__ import annotations

import csv
import gzip
import io
import subprocess
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIR = ROOT / "data" / "geo_myotube"
SIGNATURE_DIR = ROOT / "signatures" / "myotube"
MANIFEST_PATH = ROOT / "results" / "geo_myotube_signatures" / "manifest.csv"


def read_geo_table(path: Path, begin: str, end: str) -> pd.DataFrame:
    with gzip.open(path, "rt") as handle:
        lines = handle.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith(begin))
    stop = next(i for i, line in enumerate(lines) if line.startswith(end))
    return pd.read_csv(io.StringIO("".join(lines[start + 1 : stop])), sep="\t",
                       low_memory=False, dtype={"ID_REF": str, "ID": str})


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    adjusted = np.full(p_values.shape, np.nan, dtype=float)
    finite = np.isfinite(p_values)
    if finite.any():
        adjusted[finite] = multipletests(p_values[finite], method="fdr_bh")[1]
    return adjusted


def paired_statistics(values: pd.DataFrame, differences: np.ndarray) -> pd.DataFrame:
    result = stats.ttest_1samp(differences, popmean=0, axis=1, nan_policy="omit")
    score = np.asarray(result.statistic, dtype=float)
    p_value = np.asarray(result.pvalue, dtype=float)
    score[~np.isfinite(score)] = 0
    p_value[~np.isfinite(p_value)] = 1
    return pd.DataFrame(
        {"score": score, "p_value": p_value, "padj": benjamini_hochberg(p_value)},
        index=values.index,
    )


def approved_symbol_map() -> dict[str, str]:
    atlas = ad.read_h5ad(ROOT / "data" / "gtex_sn_skeletalmuscle.h5ad", backed="r")
    try:
        mapping = pd.DataFrame({
            "ensembl": atlas.var["gene_ids"].astype(str).str.split(".").str[0].to_numpy(),
            "symbol": atlas.var["Approved symbol"].astype("string").str.strip().to_numpy(),
        })
    finally:
        atlas.file.close()
    mapping = mapping.dropna().drop_duplicates()
    mapping = mapping[~mapping["symbol"].isin(["", "nan", "None"])]
    unique = mapping.groupby("ensembl")["symbol"].nunique()
    mapping = mapping[mapping["ensembl"].isin(unique[unique == 1].index)]
    return mapping.drop_duplicates("ensembl").set_index("ensembl")["symbol"].to_dict()


def export_significant(
    stats_table: pd.DataFrame,
    *,
    filename: str,
    accession: str,
    treatment: str,
    contrast: str,
    n_donors: int,
    tested_genes: int,
    score_method: str,
    source_file: str,
    source_url: str,
    cd36_p_value: float | None = None,
    cd36_padj: float | None = None,
    notes: str = "",
) -> dict:
    significant = stats_table.loc[stats_table["padj"] < 0.05].copy()
    significant = significant.loc[np.isfinite(significant["score"])]
    significant = significant.loc[significant["score"] != 0]
    significant["direction"] = np.where(significant["score"] > 0, "up", "down")
    output = SIGNATURE_DIR / filename
    if len(significant):
        significant.assign(gene=significant.index.astype(str))[
            ["gene", "score", "direction"]
        ].assign(abs_score=lambda frame: frame["score"].abs()).sort_values(
            "abs_score", ascending=False
        ).drop(columns="abs_score").to_csv(
            output, index=False, float_format="%.8g"
        )
    elif output.exists():
        output.unlink()

    return {
        "accession": accession,
        "treatment": treatment,
        "contrast": contrast,
        "n_donors": n_donors,
        "n_tested_genes": tested_genes,
        "n_fdr05": int((stats_table["padj"] < 0.05).sum()),
        "n_exported_genes": len(significant),
        "score_method": score_method,
        "signature_file": str(output.relative_to(ROOT)) if len(significant) else "",
        "cd36_score": stats_table.loc["CD36", "score"] if "CD36" in stats_table.index else np.nan,
        "cd36_p_value": cd36_p_value,
        "cd36_padj": cd36_padj,
        "source_file": source_file,
        "source_url": source_url,
        "notes": notes,
    }


def metformin_signatures(symbol_map: dict[str, str]) -> list[dict]:
    path = SOURCE_DIR / "GSE229658_Linear_normalized_data.txt.gz"
    raw = pd.read_csv(path, sep="\t", index_col=0)
    ensembl = raw.index.astype(str).str.split(".").str[0]
    raw.index = ensembl
    raw = raw.apply(pd.to_numeric, errors="raise")
    raw = raw.loc[raw.index.isin(symbol_map)]
    raw.index = raw.index.map(symbol_map)
    log_expression = np.log2(raw.clip(lower=0) + 1)
    log_expression = log_expression.groupby(level=0).mean()
    donor_ids = sorted({column.split("_", 1)[1] for column in raw.columns})
    controls = [f"M0_{donor}" for donor in donor_ids]
    if any(column not in log_expression for column in controls):
        raise ValueError("GSE229658 control samples do not match the dose-series donor IDs.")

    records = []
    for dose in ["M1", "M2", "M3", "M4", "M5"]:
        treated = [f"{dose}_{donor}" for donor in donor_ids]
        if any(column not in log_expression for column in treated):
            raise ValueError(f"GSE229658 is missing a sample for {dose}.")
        differences = (
            log_expression[treated].to_numpy() - log_expression[controls].to_numpy()
        )
        table = paired_statistics(log_expression, differences)
        dose_um = {"M1": 16, "M2": 39, "M3": 78, "M4": 388, "M5": 776}[dose]
        records.append(export_significant(
            table,
            filename=f"gse229658_metformin_{dose_um}uM_48h.csv",
            accession="GSE229658",
            treatment=f"metformin_{dose_um}uM_48h",
            contrast=f"{dose_um} uM vs 0 uM, 48 h",
            n_donors=len(donor_ids),
            tested_genes=len(table),
            score_method="paired t statistic on log2(normalized expression + 1)",
            source_file=path.name,
            source_url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE229nnn/GSE229658/suppl/GSE229658_Linear_normalized_data.txt.gz",
            cd36_p_value=table.loc["CD36", "p_value"] if "CD36" in table.index else None,
            cd36_padj=table.loc["CD36", "padj"] if "CD36" in table.index else None,
            notes="GEO sample titles say primary human myoblasts; the series summary describes human myotubes.",
        ))
    return records


def cardarine_signature() -> dict:
    matrix_path = SOURCE_DIR / "GSE40789_series_matrix.txt.gz"
    annotation_path = SOURCE_DIR / "GPL6244.annot.gz"
    with gzip.open(matrix_path, "rt") as handle:
        matrix_lines = handle.readlines()
    title_line = next(line for line in matrix_lines if line.startswith("!Sample_title\t"))
    sample_titles = next(csv.reader([title_line], delimiter="\t"))[1:]
    matrix = read_geo_table(
        matrix_path, "!series_matrix_table_begin", "!series_matrix_table_end"
    ).set_index("ID_REF")
    matrix = matrix.apply(pd.to_numeric, errors="raise")
    matrix.columns = sample_titles

    annotation = read_geo_table(
        annotation_path, "!platform_table_begin", "!platform_table_end"
    ).set_index("ID")
    symbols = annotation["Gene symbol"].astype("string").str.strip()
    single_symbol = symbols.str.split("///").str.len().eq(1)
    probe_to_symbol = symbols[single_symbol].dropna()
    matrix = matrix.loc[matrix.index.isin(probe_to_symbol.index)].copy()
    matrix.index = matrix.index.map(probe_to_symbol)
    expression = matrix.groupby(level=0).mean()

    donor_ids = sorted({title.split("_", 1)[0] for title in sample_titles})
    controls = [f"{donor}_ctrl" for donor in donor_ids]
    treated = [f"{donor}_GW501516" for donor in donor_ids]
    if any(column not in expression for column in controls + treated):
        raise ValueError("GSE40789 sample titles do not match the expected paired design.")
    with tempfile.TemporaryDirectory(prefix="gse40789_limma_") as temporary_dir:
        temporary = Path(temporary_dir)
        expression_path = temporary / "expression.tsv"
        result_path = temporary / "limma.tsv"
        expression.to_csv(expression_path, sep="\t", index_label="gene")
        subprocess.run(
            [
                "Rscript",
                str(ROOT / "notebooks" / "06_gse40789_limma.R"),
                str(expression_path),
                str(result_path),
            ],
            check=True,
        )
        table = pd.read_csv(result_path, sep="\t").set_index("gene")
        table = table.rename(columns={"t": "score", "P.Value": "p_value", "adj.P.Val": "padj"})
    cd36 = table.loc["CD36"] if "CD36" in table.index else None
    return export_significant(
        table,
        filename="gse40789_cardarine_gw501516_10nM_96h.csv",
        accession="GSE40789",
        treatment="GW501516 (Cardarine), 10 nM, 96 h",
        contrast="GW501516 vs vehicle, 96 h",
        n_donors=len(donor_ids),
        tested_genes=len(table),
        score_method="limma moderated t statistic; donor-blocked paired design",
        source_file=f"{matrix_path.name}; {annotation_path.name}",
        source_url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40789/matrix/GSE40789_series_matrix.txt.gz",
        cd36_p_value=cd36["p_value"] if cd36 is not None else None,
        cd36_padj=cd36["padj"] if cd36 is not None else None,
        notes="The GEO series is GW501516, commonly called Cardarine; only unambiguous platform probe-to-symbol mappings were used. Requires R and Bioconductor limma.",
    )


def statin_signatures() -> list[dict]:
    records = []
    definitions = [
        ("GSE107998_Simvastatin_DEG.txt.gz", "simvastatin", "Simvastatin"),
        ("GSE107998_Rosuvastatin_DEG.txt.gz", "rosuvastatin", "Rosuvastatin"),
    ]
    for source_name, treatment, label in definitions:
        path = SOURCE_DIR / source_name
        data = pd.read_csv(path, sep="\t").set_index("id")
        data = data.loc[~data.index.duplicated(keep="first")]
        log_fc = pd.to_numeric(data["log2FoldChange"], errors="coerce")
        p_value = pd.to_numeric(data["pval"], errors="coerce")
        padj = pd.to_numeric(data["padj"], errors="coerce")
        valid = np.isfinite(log_fc) & np.isfinite(p_value) & (p_value >= 0) & (p_value <= 1)
        data = data.loc[valid].copy()
        log_fc = log_fc.loc[valid]
        p_value = p_value.loc[valid]
        padj = padj.loc[valid]
        z_score = norm.isf(np.maximum(p_value.to_numpy(), 1e-300) / 2)
        z_score = np.minimum(z_score, 38)
        score = np.sign(log_fc.to_numpy()) * z_score
        table = pd.DataFrame({
            "score": score,
            "p_value": p_value.to_numpy(),
            "padj": padj.to_numpy(),
            "log2_fold_change": log_fc.to_numpy(),
        }, index=data.index.astype(str))
        cd36_p = table.loc["CD36", "p_value"] if "CD36" in table.index else None
        cd36_q = table.loc["CD36", "padj"] if "CD36" in table.index else None
        records.append(export_significant(
            table,
            filename=f"gse107998_{treatment}_5uM.csv",
            accession="GSE107998",
            treatment=f"{label}, 5 uM",
            contrast=f"{label} vs DMSO, 5 uM",
            n_donors=25,
            tested_genes=len(table),
            score_method="signed normal z from two-sided DESeq2 p value; sign from log2 fold change",
            source_file=path.name,
            source_url=f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE107nnn/GSE107998/suppl/{source_name}",
            cd36_p_value=cd36_p,
            cd36_padj=cd36_q,
            notes="Exported only genes with GEO-supplied padj < 0.05.",
        ))
    return records


def main() -> None:
    SIGNATURE_DIR.mkdir(parents=True, exist_ok=True)
    MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)
    symbols = approved_symbol_map()
    records = (
        metformin_signatures(symbols)
        + [cardarine_signature()]
        + statin_signatures()
    )
    pd.DataFrame(records).to_csv(MANIFEST_PATH, index=False)
    print(pd.DataFrame(records)[
        ["accession", "treatment", "n_tested_genes", "n_fdr05",
         "n_exported_genes", "signature_file", "cd36_padj"]
    ].to_string(index=False))
    print(f"\nManifest: {MANIFEST_PATH}")


if __name__ == "__main__":
    main()
