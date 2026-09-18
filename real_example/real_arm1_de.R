#!/usr/bin/env Rscript
# =============================================================================
# REAL-data example — Arm 1 bulk DE on a real fission-yeast RNA-seq experiment
# (S. pombe wt-vs-mutant stress timecourse; Data Carpentry / Bahler lab).
# Same statistical core as arm1_bulk/R/deseq2_edger.R (DESeq2 lfcShrink + edgeR
# TMM/glmQLF + concordance), entering at the count-matrix level (featureCounts
# path). Estimated counts are rounded to integers. Two contrasts:
#   (A) stress response: wild-type 0 vs 180 min (3 vs 3)
#   (B) strain effect:   mutant vs wild-type, all timepoints (~minute+strain, 18 vs 18)
# Usage:  Rscript real_example/real_arm1_de.R   (after fetch_real_data.sh)
# =============================================================================
suppressPackageStartupMessages({ library(DESeq2); library(edgeR); library(data.table); library(ggplot2) })

args <- commandArgs(trailingOnly = FALSE)
self <- sub("^--file=", "", args[grep("^--file=", args)])
here <- normalizePath(dirname(self))
DATA <- file.path(here, "data"); OUT <- file.path(here, "results"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cts_all <- as.data.frame(fread(file.path(DATA, "counts_raw.csv"))); rownames(cts_all) <- cts_all$gene; cts_all$gene <- NULL
info <- as.data.frame(fread(file.path(DATA, "sample_info.csv")))
info <- info[match(colnames(cts_all), info$sample), ]
info$strain <- factor(info$strain, levels = c("wt","mut")); info$minute <- factor(info$minute)

# ---- shared concordance + plotting -----------------------------------------
run_pair <- function(cts, coldata, design, coef, title, tag, lfc_cut) {
  dds <- DESeqDataSetFromMatrix(cts, colData = coldata, design = design)
  dds <- dds[rowSums(counts(dds)) >= 10, ]; dds <- DESeq(dds, quiet = TRUE)
  rd <- as.data.frame(lfcShrink(dds, coef = coef, type = "normal")); rd$gene <- rownames(rd)
  dm <- model.matrix(design, coldata)
  y <- DGEList(cts); keep <- filterByExpr(y, dm); y <- y[keep,,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y); y <- estimateDisp(y, dm)
  fit <- glmQLFit(y, dm); qlf <- glmQLFTest(fit, coef = colnames(dm)[coef])
  re <- topTags(qlf, n = Inf)$table; re$gene <- rownames(re)
  m <- merge(rd[,c("gene","log2FoldChange","padj")], re[,c("gene","logFC","FDR")], by = "gene")
  setDT(m); setnames(m, c("gene","lfc_deseq","padj_deseq","lfc_edger","padj_edger"))
  sd <- m[!is.na(padj_deseq) & padj_deseq<0.05 & abs(lfc_deseq)>lfc_cut, gene]
  se <- m[!is.na(padj_edger) & padj_edger<0.05 & abs(lfc_edger)>lfc_cut, gene]
  inter <- length(intersect(sd,se)); uni <- max(length(union(sd,se)),1); r <- cor(m$lfc_deseq, m$lfc_edger, use="complete.obs")
  m[, agree := fifelse(gene %in% intersect(sd,se),"both", fifelse(gene %in% sd,"DESeq2 only", fifelse(gene %in% se,"edgeR only","n.s.")))]
  p <- ggplot(m, aes(lfc_deseq, lfc_edger, colour=agree)) +
    geom_abline(slope=1, intercept=0, linetype=2, colour="grey50") + geom_point(size=.7, alpha=.6) +
    scale_colour_manual(values=c("both"="#E15759","DESeq2 only"="#4E79A7","edgeR only"="#59A14F","n.s."="#BAB0AC")) +
    labs(title=title, subtitle=sprintf("DESeq2=%d edgeR=%d shared=%d Jaccard=%.2f r=%.3f (padj<0.05, |lfc|>%.2g)",
         length(sd),length(se),inter,inter/uni,r,lfc_cut), x="log2FC (DESeq2)", y="log2FC (edgeR)", colour=NULL) + theme_bw(base_size=12)
  ggsave(file.path(OUT, paste0("real_", tag, "_concordance.png")), p, width=7, height=5.5, dpi=150)
  fwrite(m[order(padj_deseq)], file.path(OUT, paste0("real_", tag, "_results.tsv")), sep="\t")
  cat(sprintf(">> [%s] DESeq2=%d edgeR=%d shared=%d Jaccard=%.2f r=%.3f\n", tag, length(sd), length(se), inter, inter/uni, r))
  list(dds=dds, deseq2=length(sd), edger=length(se), shared=inter, jaccard=round(inter/uni,3), r=round(r,3))
}

# ---- (A) stress response: wt 0 vs 180 --------------------------------------
selA <- info[info$strain=="wt" & info$minute %in% c(0,180), ]
selA$condition <- factor(ifelse(selA$minute==0,"t0","t180"), levels=c("t0","t180"))
ctsA <- round(as.matrix(cts_all[, selA$sample])); mode(ctsA) <- "integer"
A <- run_pair(ctsA, data.frame(condition=selA$condition, row.names=selA$sample), ~condition, 2,
              "Real yeast: wild-type stress response (0 vs 180 min)", "stress", 1)

# ---- (B) strain effect: mut vs wt across timecourse ------------------------
cts <- round(as.matrix(cts_all)); mode(cts) <- "integer"
B <- run_pair(cts, info, ~minute + strain, which(colnames(model.matrix(~minute+strain, info))=="strainmut"),
              "Real yeast: mutant vs wild-type (all timepoints, ~minute+strain)", "strain", 0)

# ---- PCA of the whole experiment -------------------------------------------
vsd <- vst(B$dds, blind=TRUE); pd <- plotPCA(vsd, intgroup=c("strain","minute"), returnData=TRUE); pv <- round(100*attr(pd,"percentVar"))
ggsave(file.path(OUT,"real_pca.png"),
  ggplot(pd, aes(PC1, PC2, colour=strain, shape=minute)) + geom_point(size=3) +
    labs(title="PCA (real yeast, 36 samples)", x=paste0("PC1: ",pv[1],"%"), y=paste0("PC2: ",pv[2],"%")) + theme_bw(),
  width=6, height=4.5, dpi=150)

writeLines(sprintf('{"stress_0v180":{"deseq2":%d,"edger":%d,"shared":%d,"jaccard":%.3f,"r":%.3f},"strain_mut_vs_wt":{"deseq2":%d,"edger":%d,"shared":%d,"jaccard":%.3f,"r":%.3f}}',
  A$deseq2,A$edger,A$shared,A$jaccard,A$r, B$deseq2,B$edger,B$shared,B$jaccard,B$r), file.path(OUT,"real_arm1_metrics.json"))
cat(">> Arm 1 real example done. Figures + tables in real_example/results/\n")
