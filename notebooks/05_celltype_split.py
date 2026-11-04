#!/usr/bin/env python3
"""Split MoTrPAC muscle signatures using GTEx nuclei and test EPS agreement.

Run with the exmim-gpu Python interpreter from any working directory.
GEO's deposited RMA expression is used rather than reprocessing the raw CELs.
All thresholds, contrasts, seeds and top-list sizes are fixed below.
"""

from pathlib import Path
import csv
import gzip
import io
import json
import subprocess
import time
from importlib.metadata import version

import anndata as ad
import gseapy as gp
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/celltype_split"
EPS_DIR = ROOT / "data/GSE44051_invitro_stimulation"
SEED = 44051
N_NULL = 1000
TOP_N = 150
SETTINGS = {"primary": (0.5, 0.25), "sensitivity": (0.6, 0.2)}
TIMEPOINTS = [("4h", "3.5_4_hr"), ("24h", "24_hr"),
              ("early", "15_30_45_min")]
MYO_TYPES = {
    "Myocyte (sk. muscle)", "Myocyte (fast-twitch)",
    "Myocyte (slow-twitch, cytoplasmic)",
    "Myocyte (fast-twitch, cytoplasmic)", "Myocyte (NMJ-rich)",
}
CLASSES = ["myonucleus-specific", "other-cell"]
START = time.monotonic()


def progress(message):
    print(f"[{(time.monotonic() - START) / 60:.1f} min] {message}", flush=True)


def geo_table(path, begin, end):
    with gzip.open(path, "rt") as handle:
        lines = list(handle)
    start = next(i for i, line in enumerate(lines) if line.startswith(begin))
    stop = next(i for i, line in enumerate(lines) if line.startswith(end))
    table = pd.read_csv(io.StringIO("".join(lines[start + 1:stop])),
                        sep="\t", dtype={0: str}, low_memory=False)
    return table, lines[:start]


def download(name, url):
    path = EPS_DIR / name
    if not path.exists():
        progress(f"Downloading public GEO file: {name}")
        temporary = path.with_suffix(path.suffix + ".part")
        subprocess.run(["curl", "--fail", "--location", "--retry", "2",
                        "--max-time", "180", url, "-o", str(temporary)], check=True)
        temporary.replace(path)
    return path


def cell_assignment():
    progress("Reading GTEx skeletal-muscle nuclei")
    atlas = ad.read_h5ad(ROOT / "data/gtex_sn_skeletalmuscle.h5ad")
    if "tissue" in atlas.obs:
        atlas = atlas[atlas.obs["tissue"].eq("skeletalmuscle")].copy()
    labels = atlas.obs["Granular cell type"]
    if labels.isna().any():
        raise ValueError("Missing granular cell-type annotations")
    counts = labels.value_counts().rename_axis("celltype").rename("n_nuclei")
    counts.to_csv(OUT / "celltype_nucleus_counts.csv")
    progress("Granular cell-type nucleus counts:\n" + counts.to_string())
    if not MYO_TYPES.issubset(set(labels)):
        raise ValueError(f"Expected myocyte labels missing: {MYO_TYPES - set(labels)}")

    approved = atlas.var["Approved symbol"].astype(str)
    valid = ~approved.isin(["nan", "", "None"])
    symbols = pd.Series(atlas.var_names, index=atlas.var_names)
    symbols.loc[valid] = approved.loc[valid]
    rename = dict(zip(atlas.var_names, symbols))
    candidates = {}
    for symbol, previous in zip(symbols, atlas.var["Previous symbols"].astype(str)):
        for old in previous.split(","):
            old = old.strip()
            if old and old != "nan":
                candidates.setdefault(old, set()).add(symbol)
    for old, targets in candidates.items():
        if old not in set(symbols) and len(targets) == 1:
            rename[old] = next(iter(targets))
    ensembl = dict(zip(atlas.var["gene_ids"].str.split(".").str[0], symbols))

    raw_counts = sparse.csr_matrix(atlas.layers["counts"], dtype=np.float64)
    if (raw_counts.data < 0).any() or not np.allclose(raw_counts.data, np.round(raw_counts.data)):
        raise ValueError("GTEx counts layer is not nonnegative integer counts")
    codes, unique_symbols = pd.factorize(symbols, sort=False)
    collapse = sparse.csr_matrix(
        (np.ones(len(codes)), (np.arange(len(codes)), codes)),
        shape=(len(codes), len(unique_symbols)))
    progress(f"Collapsing {len(symbols) - len(unique_symbols)} duplicate canonical features")
    original_features = atlas.n_vars
    atlas = ad.AnnData(X=raw_counts @ collapse, obs=atlas.obs.copy(),
                      var=pd.DataFrame(index=pd.Index(unique_symbols, name="gene")))
    symbols = pd.Series(unique_symbols, index=unique_symbols)
    if np.any(np.asarray(atlas.X.sum(axis=1)).ravel() == 0):
        raise ValueError("GTEx has zero-library nuclei")
    sc.pp.normalize_total(atlas, target_sum=1e4)
    sc.pp.log1p(atlas)
    groups = labels.astype(str).where(~labels.isin(MYO_TYPES), "Myonuclei")
    group_names = ["Myonuclei"] + sorted(set(groups) - {"Myonuclei"})
    means, fractions = {}, {}
    for group in group_names:
        block = atlas.X[np.asarray(groups == group)]
        means[group] = np.asarray(block.mean(axis=0)).ravel()
        fractions[group] = np.asarray((block > 0).mean(axis=0)).ravel()
    mean = pd.DataFrame(means, index=pd.Index(symbols, name="gene"))
    frac = pd.DataFrame(fractions, index=mean.index)
    total = mean.sum(axis=1)
    assignment = pd.DataFrame(index=mean.index)
    assignment["myo_share"] = mean["Myonuclei"].div(total.where(total > 0))
    assignment["dominant_other_celltype"] = mean.drop(columns="Myonuclei").idxmax(axis=1)
    assignment.loc[mean.drop(columns="Myonuclei").max(axis=1) == 0,
                   "dominant_other_celltype"] = ""
    assignment["fraction_myonuclei"] = frac["Myonuclei"]
    assignment["fraction_other_max"] = frac.drop(columns="Myonuclei").max(axis=1)
    assignment["detected_snRNA"] = total > 0
    for group in group_names:
        assignment[f"mean__{group}"] = mean[group]
        assignment[f"fraction__{group}"] = frac[group]
    for setting, (myo_min, other_max) in SETTINGS.items():
        name = "class" if setting == "primary" else "class_sensitivity"
        assignment[name] = np.select(
            [(assignment.myo_share >= myo_min) & (assignment.fraction_myonuclei >= 0.1),
             (assignment.myo_share < other_max) & (assignment.fraction_other_max >= 0.1)],
            ["myonucleus-specific", "other-cell"], default="ambiguous")
    assignment.to_csv(OUT / "gene_celltype_assignment.csv")
    mean.to_csv(OUT / "celltype_mean_expression.csv")
    progress("Assignments (primary):\n" + assignment["class"].value_counts().to_string())
    metadata = {
        "n_nuclei": atlas.n_obs, "n_atlas_genes": atlas.n_vars,
        "n_original_atlas_features": original_features,
        "myonuclei": int((groups == "Myonuclei").sum()),
        "myonucleus_labels": sorted(MYO_TYPES), "n_other_celltypes": len(group_names) - 1,
        "normalization": "counts layer: normalize_total(1e4), log1p",
        "gene_symbols": "GTEx approved symbol when available; original symbol otherwise. "
                        "Unambiguous previous-symbol mappings only; no alias expansion.",
    }
    return assignment, rename, ensembl, metadata


def eps_signature(rename):
    progress("Preparing paired EPS signature from GEO-deposited RMA expression")
    matrix_path = download(
        "GSE44051_series_matrix.txt.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE44nnn/GSE44051/matrix/"
        "GSE44051_series_matrix.txt.gz")
    annotation_path = download(
        "GPL6244.annot.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6244/annot/GPL6244.annot.gz")
    expression, header = geo_table(matrix_path, "!series_matrix_table_begin",
                                   "!series_matrix_table_end")
    expression = expression.set_index("ID_REF")
    annotation, _ = geo_table(annotation_path, "!platform_table_begin",
                              "!platform_table_end")
    annotation = annotation.set_index("ID")
    metadata_rows = [next(csv.reader([line], delimiter="\t")) for line in header
                     if line.startswith("!Sample_characteristics_ch1")]
    donors = next(row[1:] for row in metadata_rows if row[1].startswith("individual id:"))
    groups = next(row[1:] for row in metadata_rows if row[1].startswith("sample group:"))
    sample_meta = pd.DataFrame({"sample": expression.columns, "donor": donors, "group": groups})
    sample_meta["stimulated"] = sample_meta.group.str.contains("stimulation", case=False)
    if len(sample_meta) != 24 or sample_meta.donor.nunique() != 12:
        raise ValueError("Expected 24 samples from 12 matched EPS/control donors")
    if not (sample_meta.groupby("donor").stimulated.agg(["sum", "count"])
            .eq([1, 2]).all().all()):
        raise ValueError("EPS donor pairing is incomplete")
    sample_meta.to_csv(OUT / "eps_sample_metadata.csv", index=False)
    if expression.min().min() < 0 or expression.max().max() > 30:
        raise ValueError("Unexpected scale for deposited log2 RMA expression")
    mapped = annotation["Gene symbol"].reindex(expression.index).fillna("").str.strip()
    unambiguous = mapped.ne("") & ~mapped.str.contains("///", regex=False)
    expression = expression.loc[unambiguous].copy()
    expression["gene"] = mapped.loc[unambiguous].map(lambda g: rename.get(g, g))
    # Collapse redundant transcript clusters before testing, without selecting on t.
    expression = expression.groupby("gene", sort=True).mean()
    pairs = []
    for _, donor in sample_meta.groupby("donor", sort=True):
        stim = donor.loc[donor.stimulated, "sample"].item()
        control = donor.loc[~donor.stimulated, "sample"].item()
        pairs.append(expression[stim] - expression[control])
    differences = np.column_stack(pairs)
    result = stats.ttest_1samp(differences, 0, axis=1)
    eps = pd.DataFrame({"gene": expression.index, "t": result.statistic,
                        "logFC": differences.mean(axis=1), "p_value": result.pvalue})
    finite = np.isfinite(eps.t)
    progress(f"EPS: {len(eps)} genes; removing {(~finite).sum()} nonfinite t values")
    eps = eps.loc[finite].copy()
    eps["adj_p_value"] = multipletests(eps.p_value, method="fdr_bh")[1]
    eps = eps.sort_values(["t", "gene"], ascending=[False, True])
    eps.to_csv(OUT / "eps_gene_ranked.csv", index=False)
    return eps.set_index("gene"), {
        "eps_method": "GEO-deposited RMA log2 expression; mean transcript-cluster "
                      "expression per unambiguously mapped gene; ordinary paired t "
                      "(stimulated minus matched control), 12 donors, df=11. "
                      "Not a reconstruction of any previously used moderated t signature.",
        "eps_unambiguous_transcript_clusters": int(unambiguous.sum()),
        "eps_finite_genes": len(eps),
        "eps_sources": [str(matrix_path.relative_to(ROOT)),
                        str(annotation_path.relative_to(ROOT))],
    }


def exercise_signatures(rename, ensembl, eps):
    progress("Selecting all six prespecified difference-in-changes contrasts")
    base = ROOT / "data/motrpac"
    de_path, = base.rglob("*t06-muscle*da_dream*")
    feature_path, = base.rglob("*t06-muscle*metadata_features*")
    de = pd.read_csv(de_path, sep="\t")
    features = pd.read_csv(feature_path, sep="\t")
    features["gene"] = features.ensembl_gene.map(ensembl).fillna(
        features.gene_symbol.map(lambda g: rename.get(g, g)))
    de = de.merge(features[["feature_id", "gene"]], on="feature_id",
                  how="left", validate="many_to_one")
    progress(f"MoTrPAC rows without symbols: {int(de.gene.isna().sum())}")
    signatures, audit = {}, []
    for label, token in TIMEPOINTS:
        for group, short in [("Endur", "endurance"), ("Resist", "resistance")]:
            contrast = (
                f"group_timepointADU{group}.post_{token} - "
                f"group_timepointADU{group}.pre_exercise - "
                f"group_timepointADUControl.post_{token} + "
                "group_timepointADUControl.pre_exercise")
            contrast_rows = de.loc[de.contrast == contrast]
            selected = contrast_rows.dropna(subset=["gene", "t"]).copy()
            if selected.empty:
                raise ValueError(f"Missing prespecified contrast: {contrast}")
            selected = selected.loc[np.isfinite(selected.t)]
            # A fixed mean across duplicate mapped genes avoids effect-based selection.
            ranked = selected.groupby("gene").t.mean().rename("motrpac_t").to_frame()
            shared = ranked.join(eps.t.rename("eps_t"), how="inner")
            shared = shared.sort_values(["motrpac_t", "gene"], ascending=[False, True])
            up = shared.loc[shared.motrpac_t > 0].head(TOP_N).index.tolist()
            down = shared.loc[shared.motrpac_t < 0].sort_values(
                ["motrpac_t", "gene"]).head(TOP_N).index.tolist()
            name = f"{short}_vs_ctrl_{label}"
            signatures[name] = {"all": shared, "up": up, "down": down}
            audit.append({"contrast": name, "original_contrast": contrast,
                          "input_features": len(contrast_rows),
                          "features_without_symbol": int(contrast_rows.gene.isna().sum()),
                          "duplicate_mapped_features": int(selected.gene.duplicated().sum()),
                          "motrpac_genes": len(ranked), "shared_genes": len(shared),
                          "top_up": len(up), "top_down": len(down)})
    pd.DataFrame(audit).to_csv(OUT / "contrast_input_audit.csv", index=False)
    return signatures


def composition(assignment, signatures):
    progress("Writing Figure 1 before agreement tests")
    rows, detection, exports = [], [], []
    export = ROOT / "signatures/exercise"
    export.mkdir(parents=True, exist_ok=True)
    for contrast, signature in signatures.items():
        all_genes = signature["all"].join(assignment, how="left")
        missing = all_genes.detected_snRNA.isna()
        all_genes["detected_snRNA"] = all_genes.detected_snRNA.eq(True)
        for column in ["class", "class_sensitivity"]:
            all_genes[column] = all_genes[column].fillna("ambiguous")
        signature["all"] = all_genes
        top = all_genes.loc[signature["up"] + signature["down"]]
        detection.append({"contrast": contrast, "top_n": len(top),
                          "top_not_in_atlas": int(missing.reindex(top.index).sum()),
                          "top_not_detected": int((~top.detected_snRNA).sum()),
                          "shared_not_detected": int((~all_genes.detected_snRNA).sum())})
        top.assign(direction=np.where(top.motrpac_t > 0, "up", "down")).to_csv(
            OUT / f"{contrast}_top_genes.csv")
        for setting in SETTINGS:
            column = "class" if setting == "primary" else "class_sensitivity"
            for direction in ["up", "down"]:
                genes = [g for g in signature[direction]
                         if all_genes.loc[g, column] == "myonucleus-specific"]
                suffix = "" if setting == "primary" else "_sensitivity"
                path = export / f"{contrast}{suffix}_myo_{direction}.txt"
                path.write_text("".join(g + "\n" for g in genes))
                exports.append({"contrast": contrast, "setting": setting,
                                "direction": direction, "n_genes": len(genes),
                                "path": str(path.relative_to(ROOT)),
                                "gsea_minimum_met": len(genes) >= 10})
            category = top[column].copy()
            other = category.eq("other-cell")
            category.loc[other] = "Other: " + top.loc[other, "dominant_other_celltype"]
            for name, count in category.value_counts().items():
                rows.append({"setting": setting, "contrast": contrast,
                             "category": name, "n": int(count),
                             "fraction": count / len(top)})
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "signature_composition.csv", index=False)
    pd.DataFrame(detection).to_csv(OUT / "snRNA_detection_audit.csv", index=False)
    pd.DataFrame(exports).to_csv(OUT / "signature_export_manifest.csv", index=False)
    categories = ["myonucleus-specific"] + sorted(
        set(table.category) - {"myonucleus-specific", "ambiguous"}) + ["ambiguous"]
    palette = dict(zip(categories, plt.get_cmap("tab20").colors))
    palette.update({"myonucleus-specific": "#276749", "ambiguous": "#cccccc"})
    fig, axes = plt.subplots(1, 2, figsize=(17, 7), sharey=True)
    for ax, setting in zip(axes, SETTINGS):
        pivot = table.loc[table.setting == setting].pivot(
            index="contrast", columns="category", values="fraction")
        pivot = pivot.reindex(index=list(signatures), columns=categories).fillna(0)
        bottom = np.zeros(len(pivot))
        for category in categories:
            ax.bar(np.arange(len(pivot)), pivot[category], bottom=bottom,
                   color=palette[category], label=category, width=0.78)
            bottom += pivot[category].to_numpy()
        a, b = SETTINGS[setting]
        ax.set_title(f"{setting.title()}: myo ≥ {a}; other < {b}")
        ax.set_xticks(np.arange(len(pivot)), [n.replace("_vs_ctrl_", "\n")
                                            for n in pivot.index], rotation=30, ha="right")
        ax.set_ylim(0, 1)
        ax.set_ylabel("Fraction of top 150 up + 150 down exercise genes")
    axes[-1].legend(loc="upper left", bbox_to_anchor=(1, 1), fontsize=8)
    fig.suptitle("GTEx nuclear cell-type assignment of MoTrPAC exercise signatures\n"
                 "Shared MoTrPAC–EPS genes; missing nuclear genes included as ambiguous")
    fig.tight_layout()
    fig.savefig(OUT / "fig1_signature_composition.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    progress("Figure 1 saved")
    return table, pd.DataFrame(detection)


def agreement(signatures, eps):
    progress("Computing full-class Spearman and 1,000 size-matched random gene sets")
    rng = np.random.default_rng(SEED)
    correlations, null_rows, gene_sets, set_audit = [], [], {}, []
    for contrast, signature in signatures.items():
        shared = signature["all"]
        for setting in SETTINGS:
            column = "class" if setting == "primary" else "class_sensitivity"
            for gene_class in CLASSES:
                selected = shared.loc[shared[column] == gene_class]
                n = len(selected)
                if n < 3:
                    raise ValueError(f"Too few genes for Spearman: {contrast}, {gene_class}, {setting}")
                rho, analytic_p = stats.spearmanr(selected.motrpac_t, selected.eps_t)
                null = np.empty(N_NULL)
                values = shared[["motrpac_t", "eps_t"]].to_numpy()
                for i in range(N_NULL):
                    draw = rng.choice(len(shared), size=n, replace=False)
                    null[i] = stats.spearmanr(values[draw, 0], values[draw, 1]).statistic
                p_upper = (1 + np.count_nonzero(null >= rho)) / (N_NULL + 1)
                p_lower = (1 + np.count_nonzero(null <= rho)) / (N_NULL + 1)
                record = {"setting": setting, "contrast": contrast, "class": gene_class,
                          "n": n, "rho": rho, "analytic_p_zero": analytic_p,
                          "n_same_sign": int((selected.motrpac_t * selected.eps_t > 0).sum()),
                          "same_sign_fraction": (selected.motrpac_t * selected.eps_t > 0).mean(),
                          "empirical_p_two_sided": min(1.0, 2 * min(p_upper, p_lower)),
                          "empirical_p_upper": p_upper, "null_median": np.median(null),
                          "null_low": np.quantile(null, .025),
                          "null_high": np.quantile(null, .975)}
                correlations.append(record)
                for i, value in enumerate(null):
                    null_rows.append({"setting": setting, "contrast": contrast,
                                      "class": gene_class, "draw": i, "rho": value})
                progress(f"{setting} {contrast} {gene_class}: n={n}, rho={rho:.4f}, "
                         f"empirical p={record['empirical_p_two_sided']:.4g}")
                for direction in ["up", "down"]:
                    genes = [g for g in signature[direction] if g in selected.index]
                    key = f"{setting}|{contrast}|{gene_class}|{direction}"
                    set_audit.append({"term": key, "setting": setting, "contrast": contrast,
                                      "class": gene_class, "direction": direction,
                                      "n_genes": len(genes),
                                      "status": "tested" if len(genes) >= 10 else "below_min_size"})
                    if len(genes) >= 10:
                        gene_sets[key] = genes
    corr = pd.DataFrame(correlations)
    corr["empirical_fdr_bh"] = multipletests(corr.empirical_p_two_sided, method="fdr_bh")[1]
    corr.to_csv(OUT / "spearman_by_class.csv", index=False)
    pd.DataFrame(null_rows).to_csv(OUT / "spearman_random_set_null.csv", index=False)
    progress(f"Running prerank GSEA: {len(gene_sets)} eligible sets, 1,000 permutations")
    if not gene_sets:
        raise ValueError("No class top-gene sets meet the fixed GSEA minimum size")
    prerank = gp.prerank(
        rnk=eps[["t"]].reset_index(), gene_sets=gene_sets,
        min_size=10, max_size=TOP_N, permutation_num=N_NULL,
        weight=1, ascending=False, threads=4, seed=SEED,
        outdir=None, no_plot=True, verbose=True)
    gsea = prerank.res2d
    audit = pd.DataFrame(set_audit)
    gsea = audit.merge(gsea, left_on="term", right_on="Term", how="left",
                       validate="one_to_one")
    if gsea.loc[gsea.status == "tested", "NES"].isna().any():
        raise ValueError("GSEA omitted an eligible prespecified set")
    gsea.to_csv(OUT / "gsea_by_class.csv", index=False)
    pd.concat([
        signature["all"].loc[signature["all"]["class"] == CLASSES[0],
                             ["motrpac_t", "eps_t", "myo_share", "class_sensitivity"]]
        .assign(contrast=contrast)
        for contrast, signature in signatures.items()
    ]).to_csv(OUT / "myonuclear_gene_statistics.csv")

    fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=True)
    positions = np.arange(len(signatures))
    for ax, setting in zip(axes, SETTINGS):
        for offset, gene_class, color in [(-.13, CLASSES[0], "#276749"),
                                           (.13, CLASSES[1], "#8055a5")]:
            block = corr.loc[(corr.setting == setting) & (corr["class"] == gene_class)]
            block = block.set_index("contrast").reindex(signatures)
            x = positions + offset
            ax.vlines(x, block.null_low, block.null_high, color=color, alpha=.4,
                      linewidth=7, label=f"{gene_class}: random-set 95% range")
            ax.scatter(x, block.null_median, color=color, marker="_", s=90)
            ax.scatter(x, block.rho, color=color, marker="o", s=45,
                       label=f"{gene_class}: observed")
            for xx, (_, row) in zip(x, block.iterrows()):
                ax.annotate(f"n={int(row.n)}", (xx, row.rho), xytext=(0, 9),
                            textcoords="offset points", ha="center", fontsize=7)
        ax.axhline(0, color="black", linewidth=.7, linestyle="--")
        ax.set_ylabel("Spearman ρ: MoTrPAC t versus EPS t")
        ax.set_title(f"{setting.title()} thresholds {SETTINGS[setting]}")
        ax.legend(fontsize=8, loc="best", ncol=2)
    axes[-1].set_xticks(positions, [n.replace("_vs_ctrl_", "\n") for n in signatures])
    fig.suptitle("Does myonuclear specificity improve EPS agreement?\n"
                 "Full gene classes; null = 1,000 equal-size sets from each shared universe")
    fig.tight_layout()
    fig.savefig(OUT / "fig2_eps_agreement_by_class.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    return corr, gsea


def write_summary(assignment, composition_table, detection, corr, gsea, metadata):
    counts = assignment["class"].value_counts()
    sensitive_counts = assignment.class_sensitivity.value_counts()
    late = corr.loc[~corr.contrast.str.endswith("_early")]
    main = late.loc[late.setting == "primary"].pivot(
        index="contrast", columns="class", values="rho")
    sensitivity = late.loc[late.setting == "sensitivity"].pivot(
        index="contrast", columns="class", values="rho")
    delta = main[CLASSES[0]] - main[CLASSES[1]]
    delta_sensitive = sensitivity[CLASSES[0]] - sensitivity[CLASSES[1]]
    direction = "higher" if delta.mean() > 0 else "lower"
    tests = gsea.loc[gsea.status == "tested"]
    significant = tests.loc[tests["FDR q-val"].astype(float) < .05]
    agreeing = significant.loc[
        ((significant.direction == "up") & (significant.NES.astype(float) > 0)) |
        ((significant.direction == "down") & (significant.NES.astype(float) < 0))]
    bullets = [
        f"GTEx: **{metadata['n_nuclei']:,} nuclei**, including "
        f"**{metadata['myonuclei']:,} pooled myonuclei**, "
        f"{metadata['n_other_celltypes']} other granular cell types, and "
        f"{metadata['n_atlas_genes']:,} measured genes. Means are equally weighted "
        "across cell types after pooling myonuclei; nuclei are not pooled across all types.",
        f"Fixed primary assignments (0.5/0.25): "
        f"**{counts.get(CLASSES[0], 0):,} myonucleus-specific**, "
        f"**{counts.get(CLASSES[1], 0):,} other-cell**, "
        f"**{counts.get('ambiguous', 0):,} ambiguous** atlas genes. "
        f"Sensitivity (0.6/0.2): {sensitive_counts.get(CLASSES[0], 0):,}, "
        f"{sensitive_counts.get(CLASSES[1], 0):,}, "
        f"{sensitive_counts.get('ambiguous', 0):,}, respectively. "
        "**Only 10 primary myonuclear genes are shared with MoTrPAC and EPS; "
        "only 3 survive the stricter setting.**",
        "All six available exercise-versus-control difference-in-changes contrasts "
        "were included (endurance/resistance × ~4 h, 24 h, early 15–45 min). "
        "Top lists are 150 positive-t and 150 negative-t shared genes, selected "
        "before splitting by class; GSEA does not refill a class with lower-ranked genes.",
    ]
    for contrast, row in main.iterrows():
        s = sensitivity.loc[contrast]
        bullets.append(
            f"**{contrast}**: full-class Spearman **myo {row[CLASSES[0]]:.3f} "
            f"vs other {row[CLASSES[1]]:.3f}** (Δ={row[CLASSES[0]] - row[CLASSES[1]]:+.3f}); "
            f"sensitivity **{s[CLASSES[0]]:.3f} vs {s[CLASSES[1]]:.3f}**.")
    bullets.append(
        f"**Conclusion:** myonucleus-specific genes have {direction} mean signed "
        f"correlation at the four late contrasts (mean Δρ={delta.mean():+.3f}; "
        f"higher in {int((delta > 0).sum())}/4; sensitivity "
        f"Δρ={delta_sensitive.mean():+.3f}, higher in "
        f"{int((delta_sensitive > 0).sum())}/4). This is a descriptive class comparison, "
        "not a formal test of the difference between correlations. "
        f"GSEA: {len(agreeing)}/{len(tests)} tested sets agree with EPS at FDR <0.05; "
        f"{len(significant) - len(agreeing)} significant sets have the opposing direction. "
        f"**All {int((gsea.status != 'tested').sum())} myonuclear sets were below "
        "the fixed 10-gene minimum, so GSEA cannot compare the classes. "
        "The correlation result is suggestive, not robustly established.**")
    text = "# Cell-type split of exercise signatures\n\n" + "\n".join(
        "- " + bullet for bullet in bullets)
    text += "\n\n## Top-signature composition and nuclear detection\n\n"
    text += "| Contrast | Myo / other / ambiguous (primary) | Not detected in snRNA / top N |\n"
    text += "|---|---:|---:|\n"
    for contrast, row in detection.set_index("contrast").iterrows():
        c = composition_table.loc[(composition_table.setting == "primary") &
                                  (composition_table.contrast == contrast)]
        myo = int(c.loc[c.category == CLASSES[0], "n"].sum())
        other = int(c.loc[c.category.str.startswith("Other:"), "n"].sum())
        ambiguous = int(c.loc[c.category == "ambiguous", "n"].sum())
        text += (f"| {contrast} | {myo} / {other} / {ambiguous} | "
                 f"{row.top_not_detected} / {row.top_n} |\n")
    text += "\n## Spearman tests (full classes, both fixed thresholds)\n\n"
    text += "| Setting | Contrast | Class | N | rho | Random-set 95% range | Empirical p | BH FDR |\n"
    text += "|---|---|---|---:|---:|---|---:|---:|\n"
    for _, r in corr.iterrows():
        text += (f"| {r.setting} | {r.contrast} | {r['class']} | {r.n} | {r.rho:.4f} | "
                 f"[{r.null_low:.4f}, {r.null_high:.4f}] | "
                 f"{r.empirical_p_two_sided:.4f} | {r.empirical_fdr_bh:.4f} |\n")
    text += "\n## Prerank GSEA on the full EPS ranked gene list\n\n"
    text += "| Setting | Contrast | Class | Direction | N | NES | FDR |\n"
    text += "|---|---|---|---|---:|---:|---:|\n"
    for _, r in gsea.iterrows():
        nes = f"{float(r.NES):.3f}" if r.status == "tested" else "not tested"
        fdr = f"{float(r['FDR q-val']):.4f}" if r.status == "tested" else "N < 10"
        text += (f"| {r.setting} | {r.contrast} | {r['class']} | {r.direction} | "
                 f"{r.n_genes} | {nes} | {fdr} |\n")
    text += """

## Methods, interpretation and limitations

**Export warning:** the primary myonuclear up/down signature files contain only
0–1 genes each; all sensitivity files are empty. They are saved faithfully in
one-gene-per-line format, but are not useful enrichment-query signatures. No
threshold was relaxed and no gene list was refilled to make them larger. See
`signature_export_manifest.csv` for every file and count. `scripts/sigcom_query.py`
is not present in this workspace, so downstream execution could not be tested.

Higher rank correlation does not necessarily mean same-direction regulation.
The `n_same_sign` and `same_sign_fraction` columns in `spearman_by_class.csv`
report that distinction; the complete 10-gene primary myonuclear statistics
are in `myonuclear_gene_statistics.csv`.

The counts layer was normalized to 10,000 counts/nucleus and log1p transformed.
Cell-type expression is the mean of log-normalized values; detection fractions use
expression >0. Myonuclei are pooled across the five labels listed in
`analysis_metadata.json`, including cytoplasmic-annotated skeletal myocytes and
NMJ-rich myocytes. Satellite cells are not myonuclei. The generic skeletal-myocyte
label cannot be resolved further into type I/II from these annotations. Other
granular types remain separate, including endothelial subtypes. Thus specificity
depends on this annotation granularity and is not evidence that a gene's exercise
response originates exclusively in that cell type.

Myo share is the pooled-myonucleus mean divided by the sum of that mean and all
other cell-type means. Specificity additionally requires >=10% myonuclear
detection; other-cell additionally requires >=10% detection in at least one other
type. Missing nuclear genes have unknown share/fractions, are explicitly marked
undetected, and remain ambiguous rather than being assigned to other-cell.
Nuclear transcript capture and the atlas's prefiltered 17,695-gene feature space
can miss genes; lack of detection is not proof of biological absence.

MoTrPAC gene IDs are mapped through its feature metadata and GTEx Ensembl IDs to
approved symbols when available. Only unambiguous previous-symbol mappings are
used to harmonize old names. Duplicate MoTrPAC mapped genes, if present, use mean
t without effect-based selection. Two duplicate canonical GTEx features (ECSCR
and TBC1D26) are summed at the count level before normalization. Unannotated
MoTrPAC features are explicitly counted in `contrast_input_audit.csv`.
Spearman uses all shared genes in each class,
not the top lists. Random sets preserve paired MoTrPAC/EPS t values, sample without
replacement from the full shared universe, and are size-matched separately per
class. The plotted interval is the 2.5–97.5 percentile random-set range, NOT an
observed-correlation confidence interval. Empirical p is two-sided, twice the
smaller tail with a +1 correction; it tests unusual agreement relative to the
shared background, NOT correlation against zero or the myo-versus-other
difference. BH adjustment covers all 24 class/contrast/threshold tests. Random sets
are not matched on expression, detection or gene–gene dependence.

GSEA uses weighted prerank (weight=1), 1,000 permutations, min_size=10, max_size=150,
seed=44051, on all finite EPS gene t values. All eligible sets from both thresholds
and all six contrasts are run together as one FDR family. Positive NES for up
sets and negative NES for down sets indicate agreement; opposite signs indicate
blunting/discordance, not exercise mimicry. Overlapping sets are not independent.
The low-resolution permutation FDR is not an exact probability; FDR=0 means below
the Monte Carlo resolution, not literally zero.

"""
    text += metadata["eps_method"] + "\n\n"
    text += ("The donor-level EPS signature is derived here from deposited processed data, "
             "not from an existing ranked table. Systemic signaling, biopsy cell mixtures, "
             "nuclear capture, myotube maturation and experimental timing can all influence "
             "agreement. This analysis does not isolate a causal explanation.\n\n")
    text += ("Reproduce: `/Users/schwepa9/miniconda3/envs/exmim-gpu/bin/python "
             "notebooks/05_celltype_split.py`.\n")
    (OUT / "summary.md").write_text(text)


def validate_results(assignment, signatures, corr, gsea):
    progress("Validating thresholds, full-universe correlations, null draws and exports")
    means = assignment.filter(regex="^mean__")
    expected_share = means["mean__Myonuclei"] / means.sum(axis=1).replace(0, np.nan)
    np.testing.assert_allclose(assignment.myo_share, expected_share, equal_nan=True)
    for setting, (myo_min, other_max) in SETTINGS.items():
        column = "class" if setting == "primary" else "class_sensitivity"
        expected = np.select([
            (assignment.myo_share >= myo_min) & (assignment.fraction_myonuclei >= .1),
            (assignment.myo_share < other_max) & (assignment.fraction_other_max >= .1),
        ], [CLASSES[0], CLASSES[1]], default="ambiguous")
        np.testing.assert_array_equal(assignment[column], expected)
    nulls = pd.read_csv(OUT / "spearman_random_set_null.csv")
    for _, row in corr.iterrows():
        column = "class" if row.setting == "primary" else "class_sensitivity"
        universe = signatures[row.contrast]["all"]
        selected = universe.loc[universe[column] == row["class"]]
        assert len(selected) == row.n
        rho = np.corrcoef(stats.rankdata(selected.motrpac_t),
                          stats.rankdata(selected.eps_t))[0, 1]
        np.testing.assert_allclose(rho, row.rho)
        null = nulls.loc[(nulls.setting == row.setting) &
                         (nulls.contrast == row.contrast) &
                         (nulls["class"] == row["class"]), "rho"]
        assert len(null) == N_NULL and np.isfinite(null).all()
        np.testing.assert_allclose(np.quantile(null, [.025, .975]),
                                   [row.null_low, row.null_high])
    manifest = pd.read_csv(OUT / "signature_export_manifest.csv")
    for _, row in manifest.iterrows():
        signature = signatures[row.contrast]
        column = "class" if row.setting == "primary" else "class_sensitivity"
        expected = [g for g in signature[row.direction]
                    if signature["all"].loc[g, column] == CLASSES[0]]
        actual = (ROOT / row.path).read_text().splitlines()
        assert actual == expected and len(actual) == row.n_genes
    for contrast, signature in signatures.items():
        universe = signature["all"]
        expected_up = universe.loc[universe.motrpac_t > 0].nlargest(TOP_N, "motrpac_t")
        expected_down = universe.loc[universe.motrpac_t < 0].nsmallest(TOP_N, "motrpac_t")
        assert set(expected_up.index) == set(signature["up"]), contrast
        assert set(expected_down.index) == set(signature["down"]), contrast
    assert len(corr) == 24 and len(gsea) == 48
    assert gsea.loc[gsea.n_genes < 10, "NES"].isna().all()
    assert gsea.loc[gsea.n_genes >= 10, "NES"].notna().all()
    assert assignment.index.is_unique
    for figure in ["fig1_signature_composition.png", "fig2_eps_agreement_by_class.png"]:
        assert (OUT / figure).stat().st_size > 10000
    progress("Validation passed")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    assignment, rename, ensembl, metadata = cell_assignment()
    eps, eps_metadata = eps_signature(rename)
    signatures = exercise_signatures(rename, ensembl, eps)
    composition_table, detection = composition(assignment, signatures)
    corr, gsea = agreement(signatures, eps)
    metadata.update(eps_metadata)
    metadata.update({"seed": SEED, "n_random_sets": N_NULL,
                     "thresholds": SETTINGS, "top_n_per_direction": TOP_N,
                     "package_versions": {name: version(name) for name in
                                          ["scanpy", "gseapy", "pandas", "scipy", "numpy"]}})
    (OUT / "analysis_metadata.json").write_text(json.dumps(metadata, indent=2))
    write_summary(assignment, composition_table, detection, corr, gsea, metadata)
    validate_results(assignment, signatures, corr, gsea)
    assert eps.index.is_unique
    assert np.allclose(composition_table.groupby(["setting", "contrast"]).fraction.sum(), 1)
    progress(f"Done. Verified outputs in {OUT}")


if __name__ == "__main__":
    main()
