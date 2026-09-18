#!/usr/bin/env Rscript
# =============================================================================
# Arm 1 — Bulk RNA-seq differential expression: DESeq2 + edgeR + concordance.
# Imports Salmon quants via tximport, runs both DE engines on the SAME counts,
# and cross-checks them (the statistics showcase).
#
# Inputs (relative to project root / CWD):
#   config/config.yaml           conditions.group_A / group_B (contrast + ref)
#   config/samples.tsv           sample sheet; uses platform == 'illumina' rows
#   data/reference/tx2gene.tsv   transcript_id <tab> gene_id [<tab> symbol]
#                                (build with scripts/make_tx2gene.sh)
#   results/arm1/salmon/<id>/quant.sf
#
# Outputs (results/arm1/de/):
#   deseq2_results.tsv, edger_results.tsv, merged_results.tsv
#   ranked_genes.rnk                       <- consumed by functional.R (GSEA)
#   pca.png, ma_deseq2.png, volcano_deseq2.png
#   deseq2_vs_edger_concordance.png
#
# Usage:  Rscript arm1_bulk/R/deseq2_edger.R
# =============================================================================

suppressPackageStartupMessages({
  library(tximport); library(DESeq2); library(edgeR)
  library(data.table); library(ggplot2)
})

# ---- tunables --------------------------------------------------------------
PADJ  <- 0.05          # significance threshold
LFC   <- 1.0           # |log2FC| threshold for calling a DEG
OUTDIR <- "results/arm1/de"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

## ---- 0. Config + sample sheet ---------------------------------------------
cfg      <- yaml::read_yaml("config/config.yaml")
grpA     <- cfg$conditions$group_A          # reference level
grpB     <- cfg$conditions$group_B          # tested-against-reference level
message(sprintf(">> Contrast: %s (test) vs %s (reference)", grpB, grpA))

samples  <- fread("config/samples.tsv")
samples  <- samples[platform == "illumina" & condition %in% c(grpA, grpB)]
stopifnot(nrow(samples) >= 4)               # need replicates per group
samples[, condition := relevel(factor(condition), ref = grpA)]

files <- file.path("results/arm1/salmon", samples$sample_id, "quant.sf")
names(files) <- samples$sample_id
if (!all(file.exists(files)))
  stop("Missing quant.sf for: ",
       paste(samples$sample_id[!file.exists(files)], collapse = ", "))

## ---- 1. tximport (transcript -> gene) -------------------------------------
t2g <- fread("data/reference/tx2gene.tsv", header = FALSE)
setnames(t2g, 1:2, c("tx", "gene"))
has_symbol <- ncol(t2g) >= 3
if (has_symbol) setnames(t2g, 3, "symbol")

# tximport's ignoreTxVersion strips versions from the quant files but NOT from
# tx2gene, so we strip the tx2gene side here too -> version-insensitive matching
# whether tx2gene came from the GENCODE FASTA or a GTF.
t2g[, tx := sub("\\.[0-9]+$", "", tx)]
txi <- tximport(files, type = "salmon",
                tx2gene = t2g[, .(tx, gene)], ignoreTxVersion = TRUE)

## ---- 2. DESeq2 -------------------------------------------------------------
coldata <- data.frame(condition = samples$condition, row.names = samples$sample_id)
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]     # pre-filter
dds <- DESeq(dds)
coef_name <- paste0("condition_", grpB, "_vs_", grpA)

# apeglm-shrunken LFC when available (preferred), else the classic normal prior.
res <- if (requireNamespace("apeglm", quietly = TRUE)) {
  lfcShrink(dds, coef = coef_name, type = "apeglm")
} else {
  message(">> apeglm not found; using type='normal' shrinkage.")
  lfcShrink(dds, coef = coef_name, type = "normal")
}
# apeglm returns the shrunken LFC WITHOUT a Wald `stat` column, so take the Wald
# statistic from the standard results() (same genes, same order) for the GSEA
# ranking. Building it explicitly also avoids a bare `stat` accidentally
# resolving to ggplot2::stat() when the column is absent.
res_std <- results(dds, name = coef_name)
stopifnot(identical(rownames(res), rownames(res_std)))
deseq_df            <- as.data.frame(res)
deseq_df$gene       <- rownames(deseq_df)
deseq_df$stat_deseq <- res_std$stat
deseq_dt <- as.data.table(deseq_df)[, .(gene,
             log2FC_deseq = log2FoldChange,
             padj_deseq   = padj,
             stat_deseq   = stat_deseq)]

## ---- 3. edgeR (canonical tximport -> edgeR offset recipe) ------------------
cts     <- txi$counts
normMat <- txi$length
normMat <- normMat / exp(rowMeans(log(normMat)))
o       <- log(calcNormFactors(cts / normMat)) + log(colSums(cts / normMat))
y       <- DGEList(cts)
y       <- scaleOffset(y, t(t(log(normMat)) + o))
keep    <- filterByExpr(y, group = samples$condition)
y       <- y[keep, ]
design  <- model.matrix(~condition, data = samples)
y       <- estimateDisp(y, design)
fit     <- glmQLFit(y, design)
qlf     <- glmQLFTest(fit, coef = 2)          # 2 = conditionB effect
edger_tt <- as.data.table(topTags(qlf, n = Inf)$table, keep.rownames = "gene")
# rename the QL F-statistic column ("F") before selecting, so we never reference
# a bare `F` (which R also reads as FALSE).
setnames(edger_tt, "F", "Fstat")
edger_dt <- edger_tt[, .(gene, log2FC_edger = logFC,
                         padj_edger = FDR, stat_edger = Fstat)]

## ---- 4. Merge + concordance ------------------------------------------------
m <- merge(deseq_dt, edger_dt, by = "gene")
if (has_symbol) m <- merge(unique(t2g[, .(gene, symbol)]), m, by = "gene", all.y = TRUE)
m[, sig_deseq := !is.na(padj_deseq) & padj_deseq < PADJ & abs(log2FC_deseq) > LFC]
m[, sig_edger := !is.na(padj_edger) & padj_edger < PADJ & abs(log2FC_edger) > LFC]

n_both <- m[, sum(sig_deseq & sig_edger)]
n_d    <- m[, sum(sig_deseq)]; n_e <- m[, sum(sig_edger)]
jacc   <- n_both / max(1, (n_d + n_e - n_both))
r_lfc  <- cor(m$log2FC_deseq, m$log2FC_edger, use = "complete.obs")
message(sprintf(">> DEGs  DESeq2=%d  edgeR=%d  shared=%d  Jaccard=%.2f  r(logFC)=%.3f",
                n_d, n_e, n_both, jacc, r_lfc))

fwrite(deseq_dt, file.path(OUTDIR, "deseq2_results.tsv"), sep = "\t")
fwrite(edger_dt, file.path(OUTDIR, "edger_results.tsv"), sep = "\t")
fwrite(m,        file.path(OUTDIR, "merged_results.tsv"), sep = "\t")

# ranked list for GSEA (DESeq2 Wald stat), one "gene<tab>stat" per line
rnk <- deseq_dt[!is.na(stat_deseq)][order(-stat_deseq), .(gene, stat_deseq)]
fwrite(rnk, file.path(OUTDIR, "ranked_genes.rnk"), sep = "\t", col.names = FALSE)

## ---- 5. Figures ------------------------------------------------------------
# (a) concordance scatter: log2FC agreement between the two engines
m[, agree := fifelse(sig_deseq & sig_edger, "both",
              fifelse(sig_deseq, "DESeq2 only",
               fifelse(sig_edger, "edgeR only", "n.s.")))]
p_conc <- ggplot(m, aes(log2FC_deseq, log2FC_edger, colour = agree)) +
  geom_abline(slope = 1, linetype = "dashed", colour = "grey60") +
  geom_point(alpha = 0.5, size = 0.9) +
  labs(title = "DESeq2 vs edgeR - log2 fold-change concordance",
       subtitle = sprintf("shared DEGs=%d  Jaccard=%.2f  r=%.3f", n_both, jacc, r_lfc),
       x = "log2FC (DESeq2, apeglm)", y = "log2FC (edgeR QL)", colour = NULL) +
  theme_bw(base_size = 12)
ggsave(file.path(OUTDIR, "deseq2_vs_edger_concordance.png"), p_conc,
       width = 6.5, height = 5.5, dpi = 150)

# (b) PCA from variance-stabilised counts
vsd <- tryCatch(vst(dds, blind = TRUE),
                error = function(e) varianceStabilizingTransformation(dds))
p_pca <- plotPCA(vsd, intgroup = "condition") + theme_bw(base_size = 12) +
  ggtitle("PCA (variance-stabilised)")
ggsave(file.path(OUTDIR, "pca.png"), p_pca, width = 6, height = 5, dpi = 150)

# (c) MA + (d) volcano (DESeq2)
png(file.path(OUTDIR, "ma_deseq2.png"), width = 1000, height = 800, res = 150)
plotMA(res, main = "MA plot (DESeq2, shrunken LFC)", ylim = c(-5, 5)); dev.off()

deseq_dt[, sig := !is.na(padj_deseq) & padj_deseq < PADJ & abs(log2FC_deseq) > LFC]
p_volc <- ggplot(deseq_dt, aes(log2FC_deseq, -log10(padj_deseq), colour = sig)) +
  geom_point(alpha = 0.5, size = 0.9) +
  scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick")) +
  geom_vline(xintercept = c(-LFC, LFC), linetype = "dashed", colour = "grey60") +
  geom_hline(yintercept = -log10(PADJ), linetype = "dashed", colour = "grey60") +
  labs(title = "Volcano (DESeq2)", x = "log2FC", y = "-log10 padj", colour = "DEG") +
  theme_bw(base_size = 12)
ggsave(file.path(OUTDIR, "volcano_deseq2.png"), p_volc, width = 6, height = 5, dpi = 150)

message(">> Arm 1 DE complete. Outputs in ", OUTDIR)
