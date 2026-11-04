#!/usr/bin/env python3
"""Draw a slide-ready overview using the saved C2S run's methods and counts."""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "cell_context_c2s"
COLORS = {
    "ink": "#183247",
    "muted": "#536777",
    "line": "#D5DFE7",
    "page": "#F4F7FA",
    "exercise": "#15785E",
    "exercise_light": "#EAF6EF",
    "drug": "#7251A3",
    "drug_light": "#F0EBF8",
    "reference": "#EDF3F9",
    "amber": "#8B5C13",
    "amber_light": "#FFF5E3",
}


def main():
    settings = json.loads((OUT / "run_settings.json").read_text())
    nuclei = pd.read_csv(OUT / "sampled_nuclei.csv")
    comparisons = pd.read_csv(OUT / "per_cell_comparisons.csv")
    coverage = pd.read_csv(OUT / "signature_coverage.csv")
    exercise_rows = coverage.loc[coverage["kind"].eq("exercise")]
    if len(exercise_rows) != 1 or comparisons["exercise"].nunique() != 1:
        raise ValueError("This figure describes the single-exercise pilot run.")
    exercise_changes = comparisons[
        ["cell_id", "exercise", "exercise_sentence_changed"]
    ].drop_duplicates()
    n_changed = int(exercise_changes["exercise_sentence_changed"].sum())
    n_nuclei = len(nuclei)
    n_comparisons = len(comparisons)
    n_c2s = int(comparisons["embedding_cosine"].notna().sum())
    n_baseline = int(comparisons["expression_weighted_correlation"].notna().sum())
    n_genes = int(exercise_rows.iloc[0]["signature_genes"])
    n_drugs = int(comparisons["drug"].nunique())

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    })
    figure = plt.figure(figsize=(16, 9), facecolor=COLORS["page"])
    axis = figure.add_axes([0, 0, 1, 1])
    axis.set(xlim=(0, 1600), ylim=(900, 0))
    axis.axis("off")

    def text(x, y, content, size=13, color="ink", weight="normal",
             ha="left", va="top", **kwargs):
        return axis.text(
            x, y, content, fontsize=size, color=COLORS.get(color, color),
            weight=weight, ha=ha, va=va, linespacing=1.35, **kwargs,
        )

    def box(x, y, width, height, face="white", edge="line", radius=14):
        patch = FancyBboxPatch(
            (x, y), width, height,
            boxstyle=f"round,pad=0,rounding_size={radius}",
            linewidth=1, facecolor=COLORS.get(face, face),
            edgecolor=COLORS.get(edge, edge),
        )
        axis.add_patch(patch)
        return patch

    def arrow(start, end, color="muted", width=1.7, style="-|>", **kwargs):
        axis.add_patch(FancyArrowPatch(
            start, end, arrowstyle=style, mutation_scale=14,
            linewidth=width, color=COLORS.get(color, color), **kwargs,
        ))

    def panel(x, width, number, title):
        box(x, 177, width, 387)
        axis.add_patch(Circle((x + 28, 207), 14, color=COLORS["ink"]))
        text(x + 28, 207, str(number), size=12, color="white",
             weight="bold", ha="center", va="center")
        text(x + 52, 197, title, size=17, weight="bold")

    text(60, 35, "Can drug signatures shift muscle-cell context in an exercise-like direction?",
         size=25, weight="bold")
    text(60, 83, "A cell-type-aware, in-silico pilot using real reference nuclei and public expression signatures",
         size=15, color="muted")
    text(60, 130, "PUBLIC DATA  →  PAIRED SIMULATED SHIFTS  →  DIRECTION + CELL-TYPE COMPARISON",
         size=11, color="muted", weight="bold")

    panel(60, 420, 1, "Use real biological inputs")
    box(80, 245, 380, 91, face="reference", edge="reference")
    text(96, 258, "GTEx skeletal-muscle nuclei", size=15, weight="bold")
    text(96, 287,
         f"{n_nuclei} nuclei · {nuclei['cell_type'].nunique()} annotated cell types\n"
         f"{settings['cells_per_type']} nuclei sampled per type",
         size=12, color="muted")
    box(80, 351, 380, 77, face="exercise_light", edge="exercise_light")
    text(96, 362, "Exercise gene signature", size=15, color="exercise", weight="bold")
    text(96, 391, f"Endurance vs control · 24 h · {n_genes} genes",
         size=11.5, color="muted")
    box(80, 443, 380, 101, face="drug_light", edge="drug_light")
    text(96, 455, f"{n_drugs} GEO drug signatures", size=15, color="drug", weight="bold")
    text(96, 484,
         "Metformin (2 doses), simvastatin,\nrosuvastatin and GW501516 / Cardarine",
         size=12, color="muted")

    panel(515, 495, 2, "Simulate paired shifts")
    box(671, 246, 184, 48, face="reference", edge="reference")
    text(763, 260, "Same reference nucleus", size=12, weight="bold", ha="center")
    branch_centers = [602, 763, 924]
    branches = [
        ("Baseline", "Unchanged", "reference", "ink"),
        ("Exercise", "Apply exercise shifts", "exercise_light", "exercise"),
        ("Drug", "Apply drug shifts", "drug_light", "drug"),
    ]
    for center, (title, subtitle, face, color) in zip(branch_centers, branches):
        arrow((763, 294), (center, 323), color=color)
        box(center - 70, 326, 140, 73, face=face, edge=face)
        text(center, 339, title, size=14, color=color, weight="bold", ha="center")
        text(center, 367, subtitle, size=9.5, color="muted", ha="center")
        arrow((center, 399), (center, 430), color=color)
    box(535, 433, 455, 63, face="reference", edge="reference")
    text(763, 444, f"Re-rank the top {settings['top_genes']} expressed genes",
         size=13, weight="bold", ha="center")
    text(763, 469, "Gene-rank sentences → C2S embeddings",
         size=12, color="muted", ha="center")
    text(763, 510, "Bounded heuristic shifts from signed gene scores;\n"
         "not measured post-treatment cells.",
         size=11, color="muted", ha="center")

    panel(1045, 495, 3, "Compare shift directions")
    text(1065, 243, "Subtract the same baseline embedding;\ncompare the two resulting shift vectors.",
         size=12, color="muted")
    origin = (1110, 362)
    axis.add_patch(Circle(origin, 5, color=COLORS["ink"]))
    arrow(origin, (1342, 299), color="exercise", width=3)
    arrow(origin, (1350, 348), color="drug", width=3)
    text(1370, 287, "Exercise shift", size=12, color="exercise", weight="bold")
    text(1370, 337, "Drug shift", size=12, color="drug", weight="bold")
    text(1110, 374, "Baseline", size=10, color="muted", ha="center")
    text(1518, 399, "Conceptual geometry, not plotted results",
         size=9, color="muted", ha="right")
    box(1065, 421, 455, 58, face="reference", edge="reference")
    text(1292, 430, r"Cosine similarity: $\Delta_{\mathrm{drug}}$ vs $\Delta_{\mathrm{exercise}}$",
         size=14, weight="bold", ha="center")
    text(1292, 457, "Average valid paired scores within each cell type",
         size=11, color="muted", ha="center")
    text(1065, 493, "Independent gene-level comparison", size=12, weight="bold")
    text(1065, 519, "Expression-weighted Pearson correlation;\n"
         "requires ≥3 shared, expressed signature genes.",
         size=11, color="muted")
    arrow((483, 374), (512, 374), color="ink", width=2)
    arrow((1013, 374), (1042, 374), color="ink", width=2)

    text(60, 591, "WHAT WE WANTED TO TEST", size=11, color="muted", weight="bold")
    goals = [
        (60, "Direction", "Do drug and exercise shifts align\nor point in opposing directions?"),
        (566, "Cell-type specificity", "Does alignment differ between\nmyocytes and other muscle-resident cells?"),
        (1072, "Added value of C2S", "Does embedding alignment agree with\nor add insight beyond gene-level similarity?"),
    ]
    for x, title, subtitle in goals:
        text(x, 620, title, size=17, weight="bold")
        text(x, 651, subtitle, size=12, color="muted")

    box(60, 727, 1480, 100, face="amber_light", edge="amber_light")
    text(80, 745, "WHAT THIS RUN\nCOULD RESOLVE", size=12, color="amber", weight="bold")
    for x, value, caption in [
        (342, f"{n_changed}/{n_nuclei}", "exercise sentences changed"),
        (645, f"{n_c2s}/{n_comparisons}", "finite C2S comparisons"),
        (925, f"{n_baseline}/{n_comparisons}", "finite gene-level comparisons"),
    ]:
        text(x, 739, value, size=24, color="amber", weight="bold")
        text(x, 785, caption, size=10.5, color="amber")
    text(1228, 750, "Limited resolution", size=14, color="amber", weight="bold")
    text(1228, 780, "Hypothesis-generating only", size=11, color="amber")
    text(60, 851,
         "Positive cosine = alignment; negative = opposition; unresolved ≠ no biology. "
         "Neither is a causal exercise-mimetic / blunting claim.",
         size=11, color="muted")
    text(60, 878,
         "Methods: C2S final-layer mean embeddings · same-nucleus shifts · signed-statistic rank proxy. "
         "Counts: saved pilot-run outputs.",
         size=9, color="muted")

    OUT.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "svg", "pdf"):
        path = OUT / f"c2s_workflow_and_test.{extension}"
        figure.savefig(path, dpi=180, facecolor=figure.get_facecolor())
        print(path)
    plt.close(figure)


if __name__ == "__main__":
    main()
