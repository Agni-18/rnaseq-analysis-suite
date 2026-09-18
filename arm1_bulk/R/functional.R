#!/usr/bin/env Rscript
# =============================================================================
# Arm 1 — Functional layer: GO/KEGG over-representation + GSEA.
# Consumes the ranked gene list written by deseq2_edger.R.
#
# Inputs:
#   results/arm1/de/ranked_genes.rnk   "ENSEMBL_gene <tab> DESeq2_stat", ranked
#   results/arm1/de/merged_results.tsv (for the DEG set used in ORA)
# Outputs (results/arm1/functional/):
#   go_bp_dotplot.png, kegg_dotplot.png            (over-representation, DEGs)
#   gsea_dotplot.png, gsea_<top>.png               (GSEA on the ranked list)
#   go_ora.tsv, kegg_ora.tsv, gsea_results.tsv
#
# Usage:  Rscript arm1_bulk/R/functional.R
# NOTE: requires clusterProfiler + org.Hs.eg.db (both in environments/env_bulk.yml).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
  library(clusterProfiler); library(org.Hs.eg.db)
})

OUTDIR <- "results/arm1/functional"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
PADJ <- 0.05; LFC <- 1.0

## ---- 0. Load DE results ----------------------------------------------------
rnk <- fread("results/arm1/de/ranked_genes.rnk", header = FALSE,
             col.names = c("gene", "stat"))
rnk[, gene := sub("\\.[0-9]+$", "", gene)]          # strip version if present
merged <- fread("results/arm1/de/merged_results.tsv")
merged[, gene := sub("\\.[0-9]+$", "", gene)]

## ---- 1. Map ENSEMBL -> ENTREZ (KEGG needs ENTREZ) --------------------------
map <- bitr(rnk$gene, fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"),
            OrgDb = org.Hs.eg.db)
setDT(map)
rnk <- merge(rnk, unique(map[, .(ENSEMBL, ENTREZID)]),
             by.x = "gene", by.y = "ENSEMBL")[!is.na(ENTREZID)]

# ranked vector for GSEA (named by ENTREZ, decreasing)
geneList <- rnk[order(-stat), setNames(stat, ENTREZID)]
geneList <- geneList[!duplicated(names(geneList))]

# DEG set for over-representation (ORA)
deg <- merged[!is.na(padj_deseq) & padj_deseq < PADJ & abs(log2FC_deseq) > LFC, gene]
deg_entrez <- unique(map[ENSEMBL %in% deg, ENTREZID])
universe   <- unique(map$ENTREZID)

save_dot <- function(obj, file, title) {
  if (is.null(obj) || nrow(as.data.frame(obj)) == 0) {
    message(">> no enrichment for: ", title); return(invisible())
  }
  ggsave(file.path(OUTDIR, file),
         dotplot(obj, showCategory = 15) + ggtitle(title),
         width = 7.5, height = 6, dpi = 150)
}

## ---- 2. Over-representation (GO BP + KEGG) ---------------------------------
ego <- enrichGO(gene = deg_entrez, universe = universe, OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID", ont = "BP",
                pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
save_dot(ego, "go_bp_dotplot.png", "GO:BP over-representation (DEGs)")
if (!is.null(ego)) fwrite(as.data.table(as.data.frame(ego)),
                          file.path(OUTDIR, "go_ora.tsv"), sep = "\t")

ekegg <- tryCatch(
  enrichKEGG(gene = deg_entrez, universe = universe,
             organism = "hsa", pAdjustMethod = "BH", qvalueCutoff = 0.05),
  error = function(e) { message(">> KEGG ORA skipped: ", conditionMessage(e)); NULL })
save_dot(ekegg, "kegg_dotplot.png", "KEGG over-representation (DEGs)")
if (!is.null(ekegg)) fwrite(as.data.table(as.data.frame(ekegg)),
                            file.path(OUTDIR, "kegg_ora.tsv"), sep = "\t")

## ---- 3. GSEA (whole ranked list) -------------------------------------------
gse <- tryCatch(
  gseGO(geneList = geneList, OrgDb = org.Hs.eg.db, ont = "BP",
        keyType = "ENTREZID", minGSSize = 10, maxGSSize = 500,
        pvalueCutoff = 0.05, verbose = FALSE),
  error = function(e) { message(">> GSEA skipped: ", conditionMessage(e)); NULL })

if (!is.null(gse) && nrow(as.data.frame(gse)) > 0) {
  ggsave(file.path(OUTDIR, "gsea_dotplot.png"),
         dotplot(gse, showCategory = 15, split = ".sign") +
           facet_grid(. ~ .sign) + ggtitle("GSEA (GO:BP)"),
         width = 9, height = 6, dpi = 150)
  fwrite(as.data.table(as.data.frame(gse)),
         file.path(OUTDIR, "gsea_results.tsv"), sep = "\t")
  # a classic running-score plot for the top term
  ggsave(file.path(OUTDIR, "gsea_top.png"),
         enrichplot::gseaplot2(gse, geneSetID = 1,
                               title = as.data.frame(gse)$Description[1]),
         width = 8, height = 5, dpi = 150)
} else {
  # ensure the Snakemake target exists even if nothing passed threshold
  ggsave(file.path(OUTDIR, "gsea_dotplot.png"),
         ggplot() + annotate("text", 0, 0, label = "No GSEA terms below FDR 0.05") +
           theme_void(), width = 6, height = 3)
}

message(">> Arm 1 functional enrichment complete. Outputs in ", OUTDIR)
