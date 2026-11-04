#!/usr/bin/env Rscript
## Step 2a. Exercise signature: top N up / top N down DEGs (adj_p_value < FDR, ranked by logFC),
## restricted to genes in the LINCS query space (landmark + best inferred). Genes renamed since
## LINCS2020 are matched by Ensembl ID (feature_id) and written under the LINCS symbol (e.g. SEPTIN7 -> SEPT7).
##
## Usage: Rscript 01_select_degs.R <deg.csv> <geneinfo_beta.txt> <outdir> <name> [fdr=0.05] [n=100]
##   deg.csv needs gene_symbol, logFC, adj_p_value; optional feature_id (Ensembl), t, p_value, AveExpr
## Writes <name>_up.txt, <name>_down.txt (LINCS query), <name>.csv (gene, score, direction, ...),
##        <name>_full_ranked.csv (all genes, ranked by t if present, else logFC)
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE)
if (length(a) < 4) stop("usage: Rscript 01_select_degs.R deg.csv geneinfo_beta.txt outdir name [fdr] [n]")
fdr <- if (length(a) > 4) as.numeric(a[5]) else 0.05
n   <- if (length(a) > 5) as.integer(a[6]) else 100

d  <- fread(a[1])
gi <- fread(a[2])[feature_space %in% c("landmark", "best inferred")]
d[, lincs_symbol := fifelse(gene_symbol %in% gi$gene_symbol, gene_symbol, NA_character_)]
if ("feature_id" %in% names(d))
  d[is.na(lincs_symbol), lincs_symbol := gi$gene_symbol[match(sub("\\..*$", "", feature_id), gi$ensembl_id)]]

s <- d[adj_p_value < fdr & !is.na(lincs_symbol)][order(adj_p_value)][!duplicated(lincs_symbol)]
sig <- rbind(head(s[logFC > 0][order(-logFC)], n)[, direction := "up"],
             head(s[logFC < 0][order(logFC)], n)[, direction := "down"])
message(sum(d$adj_p_value < fdr), " genes with FDR < ", fdr, ", ", nrow(s), " in the LINCS gene space; kept ",
        sum(sig$direction == "up"), " up / ", sum(sig$direction == "down"), " down")

dir.create(a[3], showWarnings = FALSE, recursive = TRUE)
out <- \(suffix) file.path(a[3], paste0(a[4], suffix))
for (dir in c("up", "down")) writeLines(sig[direction == dir, lincs_symbol], out(paste0("_", dir, ".txt")))
fwrite(sig[, .(gene = lincs_symbol, score = logFC, direction, original_symbol = gene_symbol, adj_p_value)], out(".csv"))
keep <- intersect(c("gene_symbol", "lincs_symbol", "logFC", "t", "p_value", "adj_p_value", "AveExpr"), names(d))
fwrite(d[order(-if ("t" %in% names(d)) t else logFC), ..keep], out("_full_ranked.csv"))
