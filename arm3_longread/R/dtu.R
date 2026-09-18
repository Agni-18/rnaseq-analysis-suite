#!/usr/bin/env Rscript
# =============================================================================
# Arm 3 — Differential Transcript Usage (DTU): DRIMSeq + two-stage stageR.
# DTU asks a different question than bulk DGE: not "is the gene up/down?" but
# "does the MIX of isoforms shift between conditions?" — the payoff of long reads.
#
# Pipeline (Love et al. 2018, F1000 RNA-seq DTU workflow):
#   DRIMSeq  dmFilter -> dmPrecision -> dmFit -> dmTest   (gene + transcript p)
#   stageR   two-stage testing: screen genes, then confirm which transcripts,
#            with overall FDR (OFDR) control on the gene set.
#
# Inputs (relative to project root / CWD):
#   config/config.yaml            conditions.group_A / group_B
#   config/samples.tsv            uses platform == 'nanopore' rows
#   results/arm3/counts/transcript_counts.tsv
#         columns: feature_id <tab> gene_id <tab> <one column per nanopore sample>
#         (run_bambu.R writes this schema; IsoQuant grouped counts reshape to it)
#
# Outputs (results/arm3/dtu/):
#   dtu_results.tsv     per-transcript: gene/tx screen+confirm adj p, significant
#   dtu_genes.tsv       significant DTU genes (passed screening stage)
#   top_switch_gene_proportions.png   isoform-usage flip for the top DTU gene
#   dtu_gene_pvalues.png              gene-level p-value histogram
#
# Usage:  Rscript arm3_longread/R/dtu.R
# =============================================================================

suppressPackageStartupMessages({
  library(DRIMSeq); library(stageR); library(data.table); library(ggplot2)
})

ALPHA  <- 0.05
OUTDIR <- "results/arm3/dtu"
CNTS   <- "results/arm3/counts/transcript_counts.tsv"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

## ---- 0. Config + samples ---------------------------------------------------
cfg  <- yaml::read_yaml("config/config.yaml")
grpA <- cfg$conditions$group_A; grpB <- cfg$conditions$group_B
message(sprintf(">> DTU contrast: %s vs %s", grpB, grpA))

samples <- fread("config/samples.tsv")
samples <- samples[platform == "nanopore" & condition %in% c(grpA, grpB)]
stopifnot(nrow(samples) >= 4)
samps <- data.frame(sample_id = samples$sample_id,
                    condition = factor(samples$condition, levels = c(grpA, grpB)))

## ---- 1. Transcript counts -> DRIMSeq object --------------------------------
cts <- fread(CNTS)
stopifnot(all(c("feature_id", "gene_id") %in% names(cts)))
stopifnot(all(samps$sample_id %in% names(cts)))
counts_df <- as.data.frame(cts[, c("feature_id", "gene_id", samps$sample_id), with = FALSE])
d <- dmDSdata(counts = counts_df, samples = samps)
message(sprintf(">> loaded %d transcripts across %d genes",
                nrow(counts_df), length(unique(counts_df$gene_id))))

# Filter: keep features/genes seen at reasonable depth in enough samples.
n <- nrow(samps); n_small <- min(table(samps$condition))
d <- dmFilter(d,
              min_samps_gene_expr    = n,       min_gene_expr    = 10,
              min_samps_feature_expr = n_small, min_feature_expr = 5,
              min_samps_feature_prop = n_small, min_feature_prop = 0.05)
message(sprintf(">> after dmFilter: %d transcripts / %d genes",
                length(counts(d)$feature_id), length(unique(counts(d)$gene_id))))

## ---- 2. DRIMSeq fit + test -------------------------------------------------
design <- model.matrix(~condition, data = DRIMSeq::samples(d))
set.seed(1)
d <- dmPrecision(d, design = design)
d <- dmFit(d, design = design)
d <- dmTest(d, coef = 2)                       # 2 = conditionB effect

res_gene <- DRIMSeq::results(d)                       # gene-level
res_txp  <- DRIMSeq::results(d, level = "feature")    # transcript-level
no.na <- function(x) ifelse(is.na(x), 1, x)           # F1000 NA handling
res_gene$pvalue <- no.na(res_gene$pvalue)
res_txp$pvalue  <- no.na(res_txp$pvalue)

## ---- 3. stageR two-stage testing (OFDR control) ----------------------------
pScreen <- res_gene$pvalue; names(pScreen) <- res_gene$gene_id
pConf   <- matrix(res_txp$pvalue, ncol = 1,
                  dimnames = list(res_txp$feature_id, "transcript"))
tx2gene <- res_txp[, c("feature_id", "gene_id")]
stageObj <- stageRTx(pScreen = pScreen, pConfirmation = pConf,
                     pScreenAdjusted = FALSE, tx2gene = tx2gene)
stageObj <- stageWiseAdjustment(stageObj, method = "dtu", alpha = ALPHA)
padj <- getAdjustedPValues(stageObj, order = FALSE, onlySignificantGenes = FALSE)
setDT(padj)  # columns: geneID, txID, gene, transcript  (gene=screen, transcript=confirm)

padj[, dtu_gene := gene < ALPHA]
padj[, dtu_transcript := transcript < ALPHA]
dtu_genes <- sort(unique(padj[dtu_gene == TRUE, geneID]))
message(sprintf(">> DTU genes (screened, OFDR<%.2f): %d | significant transcripts: %d",
                ALPHA, length(dtu_genes), padj[, sum(dtu_transcript, na.rm = TRUE)]))

## ---- 4. Outputs ------------------------------------------------------------
setnames(padj, c("geneID","txID","gene","transcript"),
               c("gene_id","feature_id","padj_gene_screen","padj_tx_confirm"))
fwrite(padj, file.path(OUTDIR, "dtu_results.tsv"), sep = "\t")
fwrite(data.table(gene_id = dtu_genes), file.path(OUTDIR, "dtu_genes.tsv"), sep = "\t")

## ---- 5. Figures ------------------------------------------------------------
# (a) isoform-usage flip for the most significant DTU gene — the signature DTU plot
if (length(dtu_genes) > 0) {
  top_gene <- padj[dtu_gene == TRUE][order(padj_gene_screen)][1, gene_id]
  p_prop <- plotProportions(d, gene_id = top_gene, group_variable = "condition") +
    ggtitle(sprintf("Isoform usage: %s (DTU gene)", top_gene))
  ggsave(file.path(OUTDIR, "top_switch_gene_proportions.png"),
         p_prop, width = 7, height = 5, dpi = 150)
}

# (b) gene-level p-value histogram
p_hist <- ggplot(res_gene, aes(pvalue)) +
  geom_histogram(boundary = 0, bins = 40, fill = "steelblue", colour = "white") +
  labs(title = "DRIMSeq gene-level p-values", x = "p-value", y = "genes") +
  theme_bw(base_size = 12)
ggsave(file.path(OUTDIR, "dtu_gene_pvalues.png"), p_hist, width = 6, height = 4.5, dpi = 150)

message(">> Arm 3 DTU complete. Outputs in ", OUTDIR)
