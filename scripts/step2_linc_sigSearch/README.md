# Step 2: exercise signature → LINCS L1000 drug matching

Turns a MoTrPAC differential-expression table into an exercise signature, scores it against every LINCS L1000
drug signature (CMap-style weighted connectivity), and summarises the result per drug.

Report: [`docs/step2_linc_sigSearch_report.html`](../../docs/step2_linc_sigSearch_report.html). Download it and open it in a browser. It
covers the DEG data, the method, the results and how to read them.

## Run

```bash
# from the repo root; R >= 4.4 with data.table, matrixStats, jsonlite and Bioconductor rhdf5
#   install.packages(c("data.table", "matrixStats", "jsonlite")); BiocManager::install("rhdf5")
bash scripts/step2_linc_sigSearch/run_step2.sh data/motrpac/ranked_24_control.csv data/lincs endurance_24h
```

This takes about 4 minutes on 1 core: the search reads 136,460 signatures from the 35 GB GCTX file. More cores help
(`--cores`, or `SLURM_CPUS_PER_TASK`).

| Script | Does | Output |
|---|---|---|
| `01_select_degs.R` | FDR < 0.05, rank by logFC, top 100 up / 100 down **among genes in the LINCS gene space** (landmark + best inferred). Renamed genes are matched by Ensembl ID. | `signatures/exercise/NAME_up.txt`, `NAME_down.txt`, `NAME.csv`, `NAME_full_ranked.csv` |
| `02_lincs_search.R` | Scores every `NAME_up/NAME_down` pair in `--sig-dir` against all LINCS drug signatures | `results/lincs/drugs.csv`, `reference_drugs.csv`, `hits.csv.gz` |
| `04_hit_targets.R` | For the top 50 mimickers and top 50 reversers: each drug's target genes and its 100 most increased / 100 most decreased LINCS genes, with the MoTrPAC logFC / t / FDR of every gene | `results/lincs/hit_targets.csv` |
| `03_report.R` | HTML report for one query (uses `report_template.html`) | `docs/step2_linc_sigSearch_report.html` |

**Data (not tracked, put under `data/`):**
- `data/lincs/`: from https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020, download
  `level5_beta_trt_cp_n720216x12328.gctx` (35 GB), `siginfo_beta.txt`, `geneinfo_beta.txt` and `compoundinfo_beta.txt`.
- The DEG table: a MoTrPAC DEA table with `gene_symbol`, `logFC`, `adj_p_value`, `t`, `p_value`, `AveExpr` and `feature_id`
  (Ensembl). The current query uses skeletal muscle, endurance vs control, 24 h, difference-in-changes.

## Method

1. **Signature.** Top 100 up / 100 down DEGs by logFC (FDR < 0.05), chosen only among genes LINCS measures or infers.
   Without that restriction about half of the top genes (mostly lncRNAs and pseudogenes) would be dropped in the search.
2. **Per LINCS signature** (one drug × cell line × dose × time, exemplars only): genes are ranked by z-score.
   A weighted Kolmogorov–Smirnov enrichment score is computed for the up set (`es_up`) and the down set (`es_down`).
3. **WTCS** = (es_up − es_down) / 2 if the two have opposite signs, else 0. Positive = mimics exercise, negative = reverses it.
4. **NCS** = WTCS divided by the mean positive (or mean negative) WTCS in that cell line.
5. **Per drug:** take the strongest dose/time per cell line, then the median NCS across cell lines. Calls are
   `mimicker` (> 0), `reverser` (< 0) and `no signal` (= 0, typically because most cell lines scored 0).

The enrichment score is computed only at the query-gene positions. This gives exactly the same value as the full
running-sum version, several times faster (about 7× in a benchmark). The scoring was checked against an independent re-implementation of the original
cmapR-based script on 3,311 real signatures; differences were at most 2e-4, from tie-breaking, with no sign changes.

## Output columns (`drugs.csv`, `reference_drugs.csv`)

| column | meaning |
|---|---|
| `median_ncs` | median NCS across cell lines; about ±1 is a typical connection, above ±1.5 is strong |
| `call` | mimicker / reverser / no signal |
| `n_cell`, `n_zero` | cell lines tested / cell lines with WTCS = 0 (up and down sets moved the same way) |
| `consistency` | share of non-zero cell lines agreeing on direction |
| `rank` | position by median NCS (1 = strongest mimicker), blank when the score is 0 (about 2/3 of drugs tie there) |
| `moa`, `target` | from `compoundinfo_beta.txt` |

`hits.csv.gz` has one row per LINCS signature: `sig_id, query, es_up, es_down, wtcs, drug, cell, dose, time, ncs`.

`hit_targets.csv` has one row per top drug × gene (about 20,000 rows). `gene_type` is `target` (annotated target; only
21 of the 100 drugs have one), `drug up` or `drug down` (the drug's 100 most increased / decreased genes in LINCS).
`lincs_z` is the drug's median LINCS z-score for the gene over its strongest signature per cell line, the signatures behind
its score; LINCS gives z-scores, not logFC. `lincs_rank` is the gene's rank within the drug's up or down list. Then come the
gene's MoTrPAC `logFC`, `t` and `adj_p_value`, summarised in `motrpac` (up / down with exercise, not significant, not tested),
plus flags `is_target` and `in_exercise_signature`.

## Main results (endurance 24 h)

- **Reversers:** half of the annotated top 50 are mTOR/PI3K inhibitors (PKI-179, GDC-0980, OSI-027, PF-04691502,
  KU-0063794, SB-2343, taselisib), plus the IGF-1R inhibitor NVP-ADW742. This is coherent: exercise activates this growth
  pathway in recovering muscle.
  The targets fit the exercise data: MTOR itself goes up with exercise (logFC +0.13, FDR 0.04) and is the target of 6 top reversers;
  IGF1R goes down (logFC −0.36, FDR 0.001) and is the target of the reverser NVP-ADW742.
- **Mimickers:** 43 of 50 are unannotated BRD compounds. The rest are broad anti-proliferative drugs (DHODH inhibitors
  brequinar and SC-12267, CDK inhibitors, daunorubicin, triptolide). Interpret with care.
- **Reference drugs:** dexamethasone is a mimicker (+0.94). Four lipophilic statins lean reverser (mevastatin −0.86 …
  cerivastatin −0.42), while the hydrophilic statins and most others score 0.

## For step 3 (myotube validation)

- **(a) Myotube drug signature vs exercise.** Use `signatures/exercise/endurance_24h_full_ranked.csv` (all genes,
  ranked by `t`; `lincs_symbol` gives the LINCS name) for Spearman correlation or prerank GSEA. Use
  `endurance_24h.csv` (`gene, score, direction`) as the up/down gene sets.
- **(b) Does LINCS find the drug again?** Compare against `results/lincs/reference_drugs.csv`, which already holds every
  step-3 drug that exists in LINCS (GW-501516, metformin, 9 statins, spironolactone, eplerenone, 2 DGAT1 inhibitors).
  Aldosterone and electrical stimulation are not in LINCS. To query LINCS with a myotube signature:
  ```bash
  Rscript scripts/step2_linc_sigSearch/01_select_degs.R myotube_metformin_DE.csv data/lincs/geneinfo_beta.txt signatures/myotube metformin_myotube
  Rscript scripts/step2_linc_sigSearch/02_lincs_search.R --sig-dir signatures/myotube --data-dir data/lincs --outdir results/lincs_myotube
  ```
  The DE table needs `gene_symbol, logFC, adj_p_value`; `feature_id` (Ensembl) is optional. Array datasets may have
  fewer than 100 genes at FDR < 0.05; pass a looser cutoff or smaller N as the 5th and 6th arguments
  (e.g. `... 0.1 50`). At least 5 up and 5 down genes must be found in LINCS.
- **Caveats.** "No signal" ≠ "no effect". LINCS uses mostly cancer cell lines at 6–24 h. The exercise signature is a
  single 24 h recovery contrast.
