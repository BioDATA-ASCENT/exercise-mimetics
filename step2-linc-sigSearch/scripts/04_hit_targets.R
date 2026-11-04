#!/usr/bin/env Rscript
## Step 2d. Genes of the top hits. For the top N mimickers and top N reversers (n_cell >= 3, consistency >= 0.75,
## ranked by median NCS) write one table with, per drug,
##   - its annotated target genes (LINCS compoundinfo), and
##   - its 100 most increased and 100 most decreased genes in LINCS (median z-score over the drug's strongest
##     signature per cell line, the same signatures used for its score),
## each with how the gene changed with exercise in MoTrPAC (logFC, t, FDR).
## Usage: Rscript 04_hit_targets.R <name> <sig_dir> <results_dir> [lincs_dir=data/lincs] [n_drugs=50] [n_genes=100]
##        -> results_dir/hit_targets.csv
suppressPackageStartupMessages({ library(data.table); library(rhdf5) })
a <- commandArgs(TRUE)
if (length(a) < 3) stop("usage: Rscript 04_hit_targets.R name sig_dir results_dir [lincs_dir] [n_drugs] [n_genes]")
NAME <- a[1]; SIG <- a[2]; RES <- a[3]; DIR <- if (length(a) > 3) a[4] else "data/lincs"
N <- if (length(a) > 4) as.integer(a[5]) else 50; NG <- if (length(a) > 5) as.integer(a[6]) else 100
GCTX <- file.path(DIR, "level5_beta_trt_cp_n720216x12328.gctx")

## Top hits and the signatures behind their scores (strongest dose/time per cell line)
st  <- fread(file.path(RES, "drugs.csv"))[query == NAME & n_cell >= 3 & consistency >= 0.75]
top <- rbind(st[order(-median_ncs)][1:N], st[order(median_ncs)][1:N])[, hit_rank := seq_len(.N), by = call]
top <- top[, .(call, hit_rank, drug, median_ncs = round(median_ncs, 3), n_cell, moa, target)]
h   <- fread(cmd = paste("zcat", shQuote(file.path(RES, "hits.csv.gz"))))[query == NAME & drug %in% top$drug]
sigs <- h[h[, .I[which.max(abs(ncs))], by = .(drug, cell)]$V1, .(drug, sig_id)]

## LINCS z-scores (genes in the query space x needed signatures)
col_id <- trimws(as.character(h5read(GCTX, "0/META/COL/id"))); row_id <- trimws(as.character(h5read(GCTX, "0/META/ROW/id")))
GI   <- fread(file.path(DIR, "geneinfo_beta.txt"))[feature_space %in% c("landmark", "best inferred")]
GI   <- unique(GI[, id := as.character(gene_id)][id %in% row_id], by = "id")
ROWS <- sort(match(GI$id, row_id)); sym <- GI$gene_symbol[match(row_id[ROWS], GI$id)]
cols <- sort(match(unique(sigs$sig_id), col_id))
Z <- h5read(GCTX, "0/DATA/0/matrix", index = list(NULL, cols))[ROWS, , drop = FALSE]; h5closeAll()
colnames(Z) <- col_id[cols]
MZ <- sapply(split(sigs$sig_id, sigs$drug), \(s) apply(Z[, s, drop = FALSE], 1, median))   # genes x drugs
rownames(MZ) <- sym

## Rows: targets + top NG up / NG down genes per drug
tg <- top[nzchar(target), .(gene = strsplit(target, ", ")[[1]], gene_type = "target"), by = drug]
pg <- rbindlist(lapply(colnames(MZ), \(d) { z <- MZ[, d]
  data.table(drug = d, gene = sym[c(head(order(-z), NG), head(order(z), NG))],
             gene_type = rep(c("drug up", "drug down"), each = NG), lincs_rank = rep(seq_len(NG), 2)) }))
out <- rbind(tg, pg, fill = TRUE)
out[, lincs_z := round(MZ[cbind(match(gene, sym), match(drug, colnames(MZ)))], 2)]   # NA if a target is not a LINCS gene
out[, is_target := mapply(\(g, t) g %in% strsplit(t, ", ")[[1]], gene, top$target[match(drug, top$drug)])]

## MoTrPAC (match on the MoTrPAC symbol or on its LINCS name, for renamed genes)
deg <- fread(file.path(SIG, paste0(NAME, "_full_ranked.csv")))
key <- unique(rbind(deg[, .(gene = gene_symbol, logFC, t, adj_p_value)],
                    deg[!is.na(lincs_symbol) & lincs_symbol != "", .(gene = lincs_symbol, logFC, t, adj_p_value)]), by = "gene")
out <- merge(out, key, by = "gene", all.x = TRUE)
sel <- fread(file.path(SIG, paste0(NAME, ".csv")))
out[, `:=`(logFC = round(logFC, 3), t = round(t, 2), adj_p_value = signif(adj_p_value, 3),
           motrpac = fifelse(is.na(logFC), "not tested", fifelse(adj_p_value < 0.05, fifelse(logFC > 0, "up with exercise", "down with exercise"), "not significant")),
           in_exercise_signature = gene %in% sel$gene)]
out <- merge(top[, !"target"], out, by = "drug")
out[, gene_type := factor(gene_type, c("target", "drug up", "drug down"))]
setcolorder(out, c("call", "hit_rank", "drug", "median_ncs", "n_cell", "moa", "gene_type", "gene", "lincs_rank", "lincs_z", "is_target",
                   "logFC", "t", "adj_p_value", "motrpac", "in_exercise_signature"))
fwrite(out[order(call, hit_rank, gene_type, lincs_rank)], file.path(RES, "hit_targets.csv"))
message("wrote hit_targets.csv: ", nrow(top), " drugs (", uniqueN(tg$drug), " with targets), ", nrow(out), " rows; ",
        round(100 * mean(out[gene_type != "target", motrpac %like% "with exercise"])), "% of drug-perturbed genes are MoTrPAC DEGs")
