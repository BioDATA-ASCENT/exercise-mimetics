#!/usr/bin/env Rscript
## Step 2b. Score query signatures against LINCS L1000 Level 5 drug signatures (CMap WTCS -> NCS)
## and summarise per drug.
##
## Usage: Rscript 02_lincs_search.R --sig-dir DIR [--data-dir data/lincs] [--outdir results/lincs] [--cores N] [--chunk 2000]
##   --sig-dir   NAME_up.txt / NAME_down.txt pairs (gene symbols, one per line); every pair is one query
##   --data-dir  LINCS2020 files: level5_beta_trt_cp_n720216x12328.gctx, siginfo_beta.txt, geneinfo_beta.txt,
##               compoundinfo_beta.txt (https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020)
## Writes (outdir): drugs.csv            one row per query x drug (>= 2 cell lines), ranked by median NCS
##                  reference_drugs.csv  positive controls and step-3 validation drugs (any number of cell lines)
##                  hits.csv.gz          one row per query x LINCS signature
## Needs R packages data.table, matrixStats, rhdf5 (Bioconductor).
suppressPackageStartupMessages({ library(data.table); library(matrixStats); library(rhdf5); library(parallel) })
a     <- commandArgs(TRUE)
opt   <- function(f, def) { i <- match(f, a); if (is.na(i)) def else a[i + 1] }
SIGS  <- opt("--sig-dir", NA); DIR <- opt("--data-dir", "data/lincs"); OUT <- opt("--outdir", "results/lincs")
CORES <- as.integer(opt("--cores", Sys.getenv("SLURM_CPUS_PER_TASK", detectCores())))
CHUNK <- as.integer(opt("--chunk", 2000))
if (is.na(SIGS)) stop("--sig-dir is required")
GCTX  <- file.path(DIR, "level5_beta_trt_cp_n720216x12328.gctx")
log_msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)
setDTthreads(1)  # parallelism comes from forked workers

## Reference drugs: positive controls and the drugs with myotube data in step 3 (LINCS cmap_name spelling)
REF <- rbind(
  data.table(group = "control: cAMP / beta-adrenergic", drug = c("forskolin", "colforsin", "rolipram", "isoprenaline", "clenbuterol")),
  data.table(group = "control: glucocorticoid", drug = "dexamethasone"),
  data.table(group = "GSE40789: PPARdelta agonist", drug = "GW-501516"),
  data.table(group = "GSE229658: metformin", drug = "metformin"),
  data.table(group = "GSE107998: lipophilic statin", drug = c("simvastatin", "atorvastatin", "lovastatin", "fluvastatin",
                                                              "cerivastatin", "pitavastatin", "mevastatin")),
  data.table(group = "GSE107998: hydrophilic statin", drug = c("pravastatin", "rosuvastatin")),
  data.table(group = "GSE70822/GSE77212: MR antagonist", drug = c("spironolactone", "eplerenone")),
  data.table(group = "GSE70822/GSE77212: DGAT1 inhibitor", drug = c("PF-04620110", "AZD-7687")))

## ---- Queries and metadata --------------------------------------------------------------------
ups <- list.files(SIGS, "_up\\.txt$", full.names = TRUE)
if (!length(ups)) stop("no *_up.txt files in ", SIGS)
rd   <- function(f) { g <- trimws(readLines(f, warn = FALSE)); unique(g[nzchar(g) & !startsWith(g, "#")]) }
sigs <- setNames(lapply(ups, \(u) list(up = rd(u), down = rd(sub("_up\\.txt$", "_down.txt", u)))),
                 sub("_up\\.txt$", "", basename(ups)))

ids    <- \(p) trimws(as.character(h5read(GCTX, p)))
col_id <- ids("0/META/COL/id"); row_id <- ids("0/META/ROW/id"); h5closeAll()
# one exemplar signature per compound x cell line x dose x time
SI <- fread(file.path(DIR, "siginfo_beta.txt"))[pert_type == "trt_cp" & is_exemplar_sig == 1 & sig_id %in% col_id,
        .(sig_id, drug = cmap_name, cell = cell_iname, dose = pert_idose, time = pert_itime, col = match(sig_id, col_id))][order(col)]
# gene universe = landmark + best inferred (the CMap query space, 10,174 genes)
GI   <- fread(file.path(DIR, "geneinfo_beta.txt"))[feature_space %in% c("landmark", "best inferred")]
GI   <- unique(GI[, id := as.character(gene_id)][id %in% row_id], by = "id")
ROWS <- sort(match(GI$id, row_id)); G <- length(ROWS)
sym  <- GI$gene_symbol[match(row_id[ROWS], GI$id)]
sets <- lapply(sigs, \(s) lapply(s, \(g) which(sym %in% g)))
log_msg(nrow(SI), " LINCS signatures, ", G, " genes")
for (q in names(sets)) {
  n <- lengths(sets[[q]]); log_msg(q, ": ", n[1], " up / ", n[2], " down genes found in LINCS")
  if (any(n < 5)) stop(q, ": need >= 5 up and >= 5 down genes in LINCS")
}
if (length(miss <- setdiff(REF$drug, SI$drug))) log_msg("reference drugs not in LINCS: ", paste(miss, collapse = ", "))
QG    <- sort(unique(unlist(sets)))                       # rows used by any query
setsq <- lapply(sets, \(s) lapply(s, match, QG))

## ---- Scoring ---------------------------------------------------------------------------------
# Weighted Kolmogorov-Smirnov enrichment score per column, evaluated only at hit positions (exactly equal
# to the full running-sum version). R: ranks of the k query genes (k x N, 1 = most up); W: their |z|.
es <- function(R, W) {
  k  <- nrow(R); o <- order(col(R), R); R <- matrix(R[o], k); W <- matrix(W[o], k)
  Ph <- colCumsums(W); Ph <- Ph / rep(Ph[k, ], each = k)     # weighted hit fraction after each hit
  Pm <- (R - seq_len(k)) / (G - k)                          # miss fraction at each hit
  mx <- colMaxs(Ph - Pm); mn <- colMins(rbind(0, Ph[-k, , drop = FALSE]) - Pm)
  ifelse(mx >= -mn, mx, mn)                                 # signed maximum deviation
}
read_cols <- function(cols) {                               # one contiguous block if dense, else indexed
  rng <- cols[1]:cols[length(cols)]
  m <- if (length(rng) <= 3 * length(cols)) h5read(GCTX, "0/DATA/0/matrix", index = list(NULL, rng))[, cols - cols[1] + 1, drop = FALSE]
       else h5read(GCTX, "0/DATA/0/matrix", index = list(NULL, cols))
  h5closeAll(); stopifnot(nrow(m) == length(row_id)); m[ROWS, , drop = FALSE]
}
score_chunk <- function(si) {
  m <- read_cols(si$col)
  R <- colRanks(-m, ties.method = "first", preserveShape = TRUE)[QG, , drop = FALSE]
  A <- abs(m[QG, , drop = FALSE])
  rbindlist(lapply(names(setsq), \(q) {
    s <- setsq[[q]]
    u <- es(R[s$up, , drop = FALSE], A[s$up, , drop = FALSE]); d <- es(R[s$down, , drop = FALSE], A[s$down, , drop = FALSE])
    data.table(sig_id = si$sig_id, query = q, es_up = u, es_down = d,
               wtcs = ifelse(sign(u) != sign(d), (u - d) / 2, 0))   # CMap WTCS: > 0 mimics, < 0 reverses
  }))
}
ch <- split(SI, ceiling(seq_len(nrow(SI)) / CHUNK))
log_msg("scoring ", length(sigs), " queries in ", length(ch), " chunks on ", CORES, " cores")
res <- mclapply(ch, score_chunk, mc.cores = CORES, mc.preschedule = FALSE)
if (any(bad <- vapply(res, inherits, NA, "try-error"))) stop(res[[which(bad)[1]]])
h <- merge(rbindlist(res), SI[, !"col"], by = "sig_id")
# NCS: WTCS scaled by the mean positive / mean negative WTCS within each query x cell line
h[, ncs := fifelse(wtcs > 0, wtcs / mean(wtcs[wtcs > 0]), fifelse(wtcs < 0, wtcs / -mean(wtcs[wtcs < 0]), 0)), by = .(query, cell)]
log_msg("scoring done")

## ---- Per-drug summary ------------------------------------------------------------------------
# strongest dose/time per cell line, then: median NCS across cell lines; n_zero = cell lines with WTCS = 0
# (up and down sets moved the same way); consistency = share of non-zero cell lines agreeing on direction
best <- h[h[, .I[which.max(abs(ncs))], by = .(query, drug, cell)]$V1]
per_drug <- best[, .(n_cell = .N, n_zero = sum(ncs == 0), median_ncs = median(ncs),
                     consistency = { nz <- ncs[ncs != 0]; if (length(nz)) max(mean(nz > 0), mean(nz < 0)) else NA_real_ },
                     cells = paste(sort(cell), collapse = ",")), by = .(query, drug)]
per_drug[, call := fifelse(median_ncs > 0, "mimicker", fifelse(median_ncs < 0, "reverser", "no signal"))]
CI <- fread(file.path(DIR, "compoundinfo_beta.txt"))[, .(
        moa    = paste(sort(unique(moa[nzchar(moa)])), collapse = "; "),
        target = paste(sort(unique(unlist(strsplit(target[nzchar(target)], "\\|")))), collapse = ", ")), by = .(drug = cmap_name)]
per_drug <- merge(per_drug, CI, by = "drug", all.x = TRUE)

drugs <- per_drug[n_cell >= 2][order(query, -median_ncs)]
drugs[, rank := fifelse(median_ncs == 0, NA_integer_, seq_len(.N)), by = query]   # ~2/3 of drugs tie at 0
ref <- merge(REF, per_drug, by = "drug")[order(query, group, -median_ncs)]

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
fwrite(drugs, file.path(OUT, "drugs.csv"))
fwrite(ref, file.path(OUT, "reference_drugs.csv"))
fwrite(h, file.path(OUT, "hits.csv.gz"))
log_msg("wrote drugs.csv, reference_drugs.csv, hits.csv.gz to ", OUT)
