#!/bin/bash
# Step 2: MoTrPAC DEG table -> exercise signature -> LINCS L1000 drug search -> HTML report.
# Run from the repo root:  bash scripts/step2_linc_sigSearch/run_step2.sh <deg.csv> [lincs_dir] [name]
# On a SLURM cluster:      sbatch -c 8 --mem 32G scripts/step2_linc_sigSearch/run_step2.sh <deg.csv> [lincs_dir] [name]
set -euo pipefail
DEG=${1:?usage: run_step2.sh deg.csv [lincs_dir=data/lincs] [name=endurance_24h]}
LINCS=${2:-data/lincs}
NAME=${3:-endurance_24h}
S=scripts/step2_linc_sigSearch

Rscript $S/01_select_degs.R "$DEG" "$LINCS/geneinfo_beta.txt" signatures/exercise "$NAME"
Rscript $S/02_lincs_search.R --sig-dir signatures/exercise --data-dir "$LINCS" --outdir results/lincs
Rscript $S/04_hit_targets.R "$NAME" signatures/exercise results/lincs "$LINCS"
Rscript $S/03_report.R "$NAME" signatures/exercise results/lincs docs/step2_linc_sigSearch_report.html
