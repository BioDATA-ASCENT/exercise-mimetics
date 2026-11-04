args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Usage: 06_gse40789_limma.R expression.tsv output.tsv")
}
if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Bioconductor package 'limma' is required for GSE40789.")
}

expression <- read.delim(args[[1]], row.names = 1, check.names = FALSE)
samples <- colnames(expression)
donor <- sub("_.*$", "", samples)
treatment <- factor(
    ifelse(grepl("_GW501516$", samples), "GW501516", "control"),
    levels = c("control", "GW501516")
)
if (any(table(donor, treatment) != 1)) {
    stop("GSE40789 must contain one control and one treated sample per donor.")
}

design <- model.matrix(~ factor(donor) + treatment)
coefficient <- match("treatmentGW501516", colnames(design))
if (is.na(coefficient)) {
    stop("Could not identify the GW501516 treatment coefficient.")
}
fit <- limma::eBayes(limma::lmFit(as.matrix(expression), design))
result <- limma::topTable(fit, coef = coefficient, number = Inf, sort.by = "none")
output <- data.frame(
    gene = rownames(result),
    t = result$t,
    P.Value = result$P.Value,
    adj.P.Val = result$adj.P.Val
)
write.table(output, args[[2]], sep = "\t", row.names = FALSE, quote = FALSE)
