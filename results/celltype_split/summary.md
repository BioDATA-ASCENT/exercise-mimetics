# Cell-type split of exercise signatures

- GTEx: **30,877 nuclei**, including **21,718 pooled myonuclei**, 17 other granular cell types, and 17,693 measured genes. Means are equally weighted across cell types after pooling myonuclei; nuclei are not pooled across all types.
- Fixed primary assignments (0.5/0.25): **13 myonucleus-specific**, **5,875 other-cell**, **11,805 ambiguous** atlas genes. Sensitivity (0.6/0.2): 3, 5,798, 11,892, respectively. **Only 10 primary myonuclear genes are shared with MoTrPAC and EPS; only 3 survive the stricter setting.**
- All six available exercise-versus-control difference-in-changes contrasts were included (endurance/resistance × ~4 h, 24 h, early 15–45 min). Top lists are 150 positive-t and 150 negative-t shared genes, selected before splitting by class; GSEA does not refill a class with lower-ranked genes.
- **endurance_vs_ctrl_24h**: full-class Spearman **myo 0.770 vs other -0.050** (Δ=+0.820); sensitivity **-0.500 vs -0.054**.
- **endurance_vs_ctrl_4h**: full-class Spearman **myo 0.927 vs other 0.094** (Δ=+0.834); sensitivity **0.500 vs 0.092**.
- **resistance_vs_ctrl_24h**: full-class Spearman **myo 0.733 vs other -0.033** (Δ=+0.767); sensitivity **1.000 vs -0.035**.
- **resistance_vs_ctrl_4h**: full-class Spearman **myo 0.891 vs other 0.053** (Δ=+0.838); sensitivity **0.500 vs 0.051**.
- **Conclusion:** myonucleus-specific genes have higher mean signed correlation at the four late contrasts (mean Δρ=+0.815; higher in 4/4; sensitivity Δρ=+0.361, higher in 3/4). This is a descriptive class comparison, not a formal test of the difference between correlations. GSEA: 6/24 tested sets agree with EPS at FDR <0.05; 4 significant sets have the opposing direction. **All 24 myonuclear sets were below the fixed 10-gene minimum, so GSEA cannot compare the classes. The correlation result is suggestive, not robustly established.**

## Top-signature composition and nuclear detection

| Contrast | Myo / other / ambiguous (primary) | Not detected in snRNA / top N |
|---|---:|---:|
| endurance_vs_ctrl_4h | 0 / 152 / 148 | 9 / 300 |
| resistance_vs_ctrl_4h | 1 / 164 / 135 | 11 / 300 |
| endurance_vs_ctrl_24h | 1 / 130 / 169 | 15 / 300 |
| resistance_vs_ctrl_24h | 0 / 155 / 145 | 8 / 300 |
| endurance_vs_ctrl_early | 0 / 119 / 181 | 25 / 300 |
| resistance_vs_ctrl_early | 0 / 122 / 178 | 20 / 300 |

## Spearman tests (full classes, both fixed thresholds)

| Setting | Contrast | Class | N | rho | Random-set 95% range | Empirical p | BH FDR |
|---|---|---|---:|---:|---|---:|---:|
| primary | endurance_vs_ctrl_4h | myonucleus-specific | 10 | 0.9273 | [-0.5758, 0.7212] | 0.0020 | 0.0480 |
| primary | endurance_vs_ctrl_4h | other-cell | 5294 | 0.0936 | [0.0661, 0.1064] | 0.4935 | 0.8152 |
| sensitivity | endurance_vs_ctrl_4h | myonucleus-specific | 3 | 0.5000 | [-1.0000, 1.0000] | 1.0000 | 1.0000 |
| sensitivity | endurance_vs_ctrl_4h | other-cell | 5221 | 0.0925 | [0.0657, 0.1062] | 0.5674 | 0.8511 |
| primary | resistance_vs_ctrl_4h | myonucleus-specific | 10 | 0.8909 | [-0.5882, 0.7091] | 0.0040 | 0.0480 |
| primary | resistance_vs_ctrl_4h | other-cell | 5294 | 0.0531 | [0.0240, 0.0641] | 0.3896 | 0.8152 |
| sensitivity | resistance_vs_ctrl_4h | myonucleus-specific | 3 | 0.5000 | [-1.0000, 1.0000] | 1.0000 | 1.0000 |
| sensitivity | resistance_vs_ctrl_4h | other-cell | 5221 | 0.0509 | [0.0221, 0.0637] | 0.4835 | 0.8152 |
| primary | endurance_vs_ctrl_24h | myonucleus-specific | 10 | 0.7697 | [-0.6364, 0.5879] | 0.0060 | 0.0480 |
| primary | endurance_vs_ctrl_24h | other-cell | 5294 | -0.0504 | [-0.0733, -0.0315] | 0.9111 | 1.0000 |
| sensitivity | endurance_vs_ctrl_24h | myonucleus-specific | 3 | -0.5000 | [-1.0000, 1.0000] | 1.0000 | 1.0000 |
| sensitivity | endurance_vs_ctrl_24h | other-cell | 5221 | -0.0539 | [-0.0721, -0.0310] | 0.8252 | 1.0000 |
| primary | resistance_vs_ctrl_24h | myonucleus-specific | 10 | 0.7333 | [-0.6245, 0.6606] | 0.0220 | 0.1319 |
| primary | resistance_vs_ctrl_24h | other-cell | 5294 | -0.0333 | [-0.0335, 0.0075] | 0.0539 | 0.2238 |
| sensitivity | resistance_vs_ctrl_24h | myonucleus-specific | 3 | 1.0000 | [-1.0000, 1.0000] | 0.3237 | 0.8152 |
| sensitivity | resistance_vs_ctrl_24h | other-cell | 5221 | -0.0349 | [-0.0354, 0.0065] | 0.0559 | 0.2238 |
| primary | endurance_vs_ctrl_early | myonucleus-specific | 10 | 0.2364 | [-0.5758, 0.6606] | 0.6274 | 0.8857 |
| primary | endurance_vs_ctrl_early | other-cell | 5294 | 0.0605 | [0.0394, 0.0816] | 0.9650 | 1.0000 |
| sensitivity | endurance_vs_ctrl_early | myonucleus-specific | 3 | 1.0000 | [-1.0000, 1.0000] | 0.3836 | 0.8152 |
| sensitivity | endurance_vs_ctrl_early | other-cell | 5221 | 0.0602 | [0.0427, 0.0832] | 0.8871 | 1.0000 |
| primary | resistance_vs_ctrl_early | myonucleus-specific | 10 | 0.3576 | [-0.6367, 0.6488] | 0.3736 | 0.8152 |
| primary | resistance_vs_ctrl_early | other-cell | 5294 | 0.0447 | [0.0310, 0.0704] | 0.5095 | 0.8152 |
| sensitivity | resistance_vs_ctrl_early | myonucleus-specific | 3 | 1.0000 | [-1.0000, 1.0000] | 0.3776 | 0.8152 |
| sensitivity | resistance_vs_ctrl_early | other-cell | 5221 | 0.0439 | [0.0310, 0.0731] | 0.4156 | 0.8152 |

## Prerank GSEA on the full EPS ranked gene list

| Setting | Contrast | Class | Direction | N | NES | FDR |
|---|---|---|---|---:|---:|---:|
| primary | endurance_vs_ctrl_4h | myonucleus-specific | up | 0 | not tested | N < 10 |
| primary | endurance_vs_ctrl_4h | myonucleus-specific | down | 0 | not tested | N < 10 |
| primary | endurance_vs_ctrl_4h | other-cell | up | 86 | 2.228 | 0.0000 |
| primary | endurance_vs_ctrl_4h | other-cell | down | 66 | 1.511 | 0.0223 |
| sensitivity | endurance_vs_ctrl_4h | myonucleus-specific | up | 0 | not tested | N < 10 |
| sensitivity | endurance_vs_ctrl_4h | myonucleus-specific | down | 0 | not tested | N < 10 |
| sensitivity | endurance_vs_ctrl_4h | other-cell | up | 82 | 2.101 | 0.0000 |
| sensitivity | endurance_vs_ctrl_4h | other-cell | down | 65 | 1.521 | 0.0217 |
| primary | resistance_vs_ctrl_4h | myonucleus-specific | up | 1 | not tested | N < 10 |
| primary | resistance_vs_ctrl_4h | myonucleus-specific | down | 0 | not tested | N < 10 |
| primary | resistance_vs_ctrl_4h | other-cell | up | 95 | 2.211 | 0.0000 |
| primary | resistance_vs_ctrl_4h | other-cell | down | 69 | 1.067 | 0.3471 |
| sensitivity | resistance_vs_ctrl_4h | myonucleus-specific | up | 0 | not tested | N < 10 |
| sensitivity | resistance_vs_ctrl_4h | myonucleus-specific | down | 0 | not tested | N < 10 |
| sensitivity | resistance_vs_ctrl_4h | other-cell | up | 89 | 2.192 | 0.0000 |
| sensitivity | resistance_vs_ctrl_4h | other-cell | down | 68 | 1.067 | 0.3260 |
| primary | endurance_vs_ctrl_24h | myonucleus-specific | up | 1 | not tested | N < 10 |
| primary | endurance_vs_ctrl_24h | myonucleus-specific | down | 0 | not tested | N < 10 |
| primary | endurance_vs_ctrl_24h | other-cell | up | 42 | 1.329 | 0.0955 |
| primary | endurance_vs_ctrl_24h | other-cell | down | 88 | 1.115 | 0.3088 |
| sensitivity | endurance_vs_ctrl_24h | myonucleus-specific | up | 0 | not tested | N < 10 |
| sensitivity | endurance_vs_ctrl_24h | myonucleus-specific | down | 0 | not tested | N < 10 |
| sensitivity | endurance_vs_ctrl_24h | other-cell | up | 41 | 1.096 | 0.3165 |
| sensitivity | endurance_vs_ctrl_24h | other-cell | down | 84 | 1.260 | 0.1424 |
| primary | resistance_vs_ctrl_24h | myonucleus-specific | up | 0 | not tested | N < 10 |
| primary | resistance_vs_ctrl_24h | myonucleus-specific | down | 0 | not tested | N < 10 |
| primary | resistance_vs_ctrl_24h | other-cell | up | 65 | -1.152 | 1.0000 |
| primary | resistance_vs_ctrl_24h | other-cell | down | 90 | 1.794 | 0.0032 |
| sensitivity | resistance_vs_ctrl_24h | myonucleus-specific | up | 0 | not tested | N < 10 |
| sensitivity | resistance_vs_ctrl_24h | myonucleus-specific | down | 0 | not tested | N < 10 |
| sensitivity | resistance_vs_ctrl_24h | other-cell | up | 62 | -1.091 | 0.8263 |
| sensitivity | resistance_vs_ctrl_24h | other-cell | down | 86 | 1.671 | 0.0063 |
| primary | endurance_vs_ctrl_early | myonucleus-specific | up | 0 | not tested | N < 10 |
| primary | endurance_vs_ctrl_early | myonucleus-specific | down | 0 | not tested | N < 10 |
| primary | endurance_vs_ctrl_early | other-cell | up | 72 | -1.061 | 0.4763 |
| primary | endurance_vs_ctrl_early | other-cell | down | 47 | -1.086 | 0.6329 |
| sensitivity | endurance_vs_ctrl_early | myonucleus-specific | up | 0 | not tested | N < 10 |
| sensitivity | endurance_vs_ctrl_early | myonucleus-specific | down | 0 | not tested | N < 10 |
| sensitivity | endurance_vs_ctrl_early | other-cell | up | 72 | -1.061 | 0.4763 |
| sensitivity | endurance_vs_ctrl_early | other-cell | down | 46 | -1.111 | 1.0000 |
| primary | resistance_vs_ctrl_early | myonucleus-specific | up | 0 | not tested | N < 10 |
| primary | resistance_vs_ctrl_early | myonucleus-specific | down | 0 | not tested | N < 10 |
| primary | resistance_vs_ctrl_early | other-cell | up | 83 | 1.780 | 0.0026 |
| primary | resistance_vs_ctrl_early | other-cell | down | 39 | -0.909 | 0.7161 |
| sensitivity | resistance_vs_ctrl_early | myonucleus-specific | up | 0 | not tested | N < 10 |
| sensitivity | resistance_vs_ctrl_early | myonucleus-specific | down | 0 | not tested | N < 10 |
| sensitivity | resistance_vs_ctrl_early | other-cell | up | 83 | 1.780 | 0.0026 |
| sensitivity | resistance_vs_ctrl_early | other-cell | down | 37 | -0.698 | 0.9482 |


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

GEO-deposited RMA log2 expression; mean transcript-cluster expression per unambiguously mapped gene; ordinary paired t (stimulated minus matched control), 12 donors, df=11. Not a reconstruction of any previously used moderated t signature.

The donor-level EPS signature is derived here from deposited processed data, not from an existing ranked table. Systemic signaling, biopsy cell mixtures, nuclear capture, myotube maturation and experimental timing can all influence agreement. This analysis does not isolate a causal explanation.

Reproduce: `/Users/schwepa9/miniconda3/envs/exmim-gpu/bin/python notebooks/05_celltype_split.py`.
