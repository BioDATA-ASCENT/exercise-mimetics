#!/usr/bin/env Rscript
## Step 2c. HTML report for one query: DEG overview (volcano, MA, p-values) + LINCS results.
## Usage: Rscript 03_report.R <name> <sig_dir> <results_dir> <out.html>
##   e.g. Rscript 03_report.R endurance_24h signatures/exercise results/lincs docs/step2_linc_sigSearch_report.html
## Needs data.table, jsonlite; uses report_template.html next to this script.
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
a <- commandArgs(TRUE)
if (length(a) < 4) stop("usage: Rscript 03_report.R name sig_dir results_dir out.html")
NAME <- a[1]; SIG <- a[2]; RES <- a[3]
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))

## DEG side
f   <- fread(file.path(SIG, paste0(NAME, "_full_ranked.csv")))
sel <- fread(file.path(SIG, paste0(NAME, ".csv")))
f[, `:=`(sig = adj_p_value < 0.05, inL = !is.na(lincs_symbol) & lincs_symbol != "")]
f[, s := sig & lincs_symbol %in% sel$gene]
deg <- f[, .(sym = gene_symbol, lfc = round(logFC, 3), nlp = round(-log10(adj_p_value), 3), ave = round(AveExpr, 2),
             sig = as.integer(sig), inL = as.integer(inL), s = as.integer(s))]
top <- rbind(f[sig == TRUE][order(-logFC)][1:15][, dir := "up"], f[sig == TRUE][order(logFC)][1:15][, dir := "down"])[
  , .(dir, gene = gene_symbol, logFC = round(logFC, 2), fdr = signif(adj_p_value, 2), inLINCS = inL)]
counts <- list(n = nrow(f), sig = sum(f$sig), up = sum(f$sig & f$logFC > 0), down = sum(f$sig & f$logFC < 0),
               sigInL = sum(f$sig & f$inL), renamed = sum(f$sig & f$inL & f$lincs_symbol != f$gene_symbol),
               cut_up = min(sel[direction == "up", score]), cut_down = max(sel[direction == "down", score]))
phist <- as.integer(table(cut(f$p_value, seq(0, 1, 0.05), include.lowest = TRUE)))

## LINCS side
h   <- fread(cmd = paste("zcat", shQuote(file.path(RES, "hits.csv.gz"))))[query == NAME]
dr  <- fread(file.path(RES, "drugs.csv"))[query == NAME]
ref <- fread(file.path(RES, "reference_drugs.csv"))[query == NAME]
strong <- dr[n_cell >= 3 & consistency >= 0.75]
cols <- \(z) z[, .(drug, n_cell, n_zero, ncs = round(median_ncs, 2), cons = round(consistency, 2), moa)]
moa_top <- function(z) { t <- sort(table(unlist(strsplit(z$moa[nzchar(z$moa)], "; "))), decreasing = TRUE)
                         head(data.table(moa = names(t), n = as.integer(t)), 8) }
br  <- seq(-2.5, 2.5, 0.1)
run <- list(n_sig = nrow(h), n_drugs = nrow(dr), n_strong = nrow(strong), frac_zero = round(mean(h$wtcs == 0), 3),
            n_nosignal = sum(dr$call == "no signal"),
            mim = cols(strong[order(-median_ncs)][1:50]), rev = cols(strong[order(median_ncs)][1:50]),
            moa_mim = moa_top(strong[order(-median_ncs)][1:50]), moa_rev = moa_top(strong[order(median_ncs)][1:50]),
            hist = as.integer(table(cut(pmin(pmax(dr$median_ncs, -2.49), 2.49), br))), hist_br = br,
            ref = merge(ref[, .(group, drug, ncs = round(median_ncs, 2), n_cell, n_zero, call)],
                        dr[, .(drug, rank)], by = "drug", all.x = TRUE)[order(group, -ncs)])

## Target genes of the top hits (from 04_hit_targets.R, optional)
targets <- NULL
if (file.exists(file.path(RES, "hit_targets.csv"))) {
  tt <- fread(file.path(RES, "hit_targets.csv"))[gene_type == "target"]
  targets <- tt[order(call, hit_rank), .(call, hit_rank, drug, ncs = round(median_ncs, 2), moa, g = gene, lfc = logFC, fdr = adj_p_value, m = motrpac)]
}

## Fill the template
data <- toJSON(list(name = NAME, counts = counts, deg = as.list(deg), phist = phist, topgenes = top, run = run, targets = targets),
               auto_unbox = TRUE, digits = NA, na = "null", dataframe = "rows")
tpl <- paste(readLines(file.path(here, "report_template.html"), warn = FALSE), collapse = "\n")
dir.create(dirname(a[4]), showWarnings = FALSE, recursive = TRUE)
writeLines(sub("/*DATA*/null", data, tpl, fixed = TRUE), a[4])
message("wrote ", a[4])
