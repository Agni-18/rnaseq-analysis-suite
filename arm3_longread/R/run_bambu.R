#!/usr/bin/env Rscript
# =============================================================================
# Arm 3 — Bambu: reference-guided transcript discovery + quantification from
# long-read (minimap2 genome) BAMs. Emits the transcript-count contract that
# dtu.R consumes:  feature_id <tab> gene_id <tab> <one column per sample>.
#
# Inputs (relative to project root):
#   data/reference/genome.subset.fa         (or genome.fa)
#   data/reference/annotation.subset.gtf    (or annotation.gtf)
#   results/arm3/bam/<sample>.bam           (from rules/align.smk minimap2)
#   config/samples.tsv                      (platform == 'nanopore' rows)
#
# Outputs:
#   results/arm3/bambu/                      native bambu output (extended GTF etc.)
#   results/arm3/counts/transcript_counts.tsv   <- DTU input contract
#
# Usage:  Rscript arm3_longread/R/run_bambu.R
# =============================================================================
suppressPackageStartupMessages({ library(bambu); library(data.table) })
`%||%` <- function(a, b) if (is.null(a)) b else a

cfg  <- yaml::read_yaml("config/config.yaml")
fa   <- if (file.exists("data/reference/genome.subset.fa")) "data/reference/genome.subset.fa" else cfg$reference$genome_fa
gtf  <- if (file.exists("data/reference/annotation.subset.gtf")) "data/reference/annotation.subset.gtf" else cfg$reference$gtf

samples <- fread("config/samples.tsv")[platform == "nanopore"]
bams <- file.path("results/arm3/bam", paste0(samples$sample_id, ".bam"))
stopifnot(all(file.exists(bams)))

annotations <- prepareAnnotations(gtf)
se <- bambu(reads = bams, annotations = annotations, genome = fa,
            ncore = cfg$threads %||% 4)
dir.create("results/arm3/bambu",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/arm3/counts", recursive = TRUE, showWarnings = FALSE)
writeBambuOutput(se, path = "results/arm3/bambu")

# ---- transcript-count contract for DTU -------------------------------------
counts <- as.matrix(assays(se)$counts)
colnames(counts) <- sub("\\.bam$", "", basename(bams))     # -> sample_id
tab <- data.table(feature_id = rownames(se),
                  gene_id    = as.character(rowData(se)$GENEID))
tab <- cbind(tab, as.data.table(round(counts)))            # bambu counts are fractional
fwrite(tab, "results/arm3/counts/transcript_counts.tsv", sep = "\t")
message(sprintf(">> wrote %d transcripts x %d samples to results/arm3/counts/transcript_counts.tsv",
                nrow(tab), ncol(counts)))

