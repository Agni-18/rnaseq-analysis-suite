#!/usr/bin/env Rscript
# =============================================================================
# Arm 1 validation: simulate a Salmon dataset with known DE ground truth, run
# the SHIPPED arm1_bulk/R/deseq2_edger.R against it, and score precision/recall/
# sign. Exits non-zero if thresholds aren't met (CI-friendly).
# Usage:  Rscript tests/validate_de.R
# =============================================================================
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = FALSE)
self <- sub("^--file=", "", args[grep("^--file=", args)])
repo <- normalizePath(file.path(dirname(self), ".."))
DE   <- file.path(repo, "arm1_bulk", "R", "deseq2_edger.R")
work <- file.path(repo, "tests", ".work", "arm1")
unlink(work, recursive = TRUE); dir.create(file.path(work, "config"), recursive = TRUE)
setwd(work)

# ---- simulate (matched to the documented Arm 1 validation) ------------------
set.seed(42)
G <- 1500L; N_DE <- 200L
samples <- c(paste0("SGNex_A549_Illumina_replicate", c(1,3,5), "_run1"),
             paste0("SGNex_K562_Illumina_replicate", c(1,3,5), "_run1"))
cond <- c("A549","A549","A549","K562","K562","K562")
gid <- sprintf("ENSG%011d.%d", 1:G, sample(1:9, G, TRUE))
ntx <- sample(1:3, G, TRUE)
tx  <- unlist(lapply(seq_len(G), function(i) sprintf("ENST%011d.%d", (i*10):(i*10+ntx[i]-1), sample(1:9, ntx[i], TRUE))))
tx_gene <- rep(gid, ntx); txlen <- sample(500:6000, length(tx), TRUE)
mu_gene <- exp(rnorm(G, 5.6, 1.6))
de_idx <- sample(G, N_DE); lfc_true <- numeric(G); lfc_true[de_idx] <- sample(c(-1,1),N_DE,TRUE)*runif(N_DE,0.6,3.2)
fold_K562 <- 2^lfc_true
prop <- unlist(lapply(ntx, function(k){p<-runif(k,.2,1); p/sum(p)}))

dir.create("data/reference", recursive = TRUE, showWarnings = FALSE)
write.table(data.frame(tx=tx, gene=tx_gene, symbol=paste0("SYM", match(tx_gene,gid))),
            "data/reference/tx2gene.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
for (s in seq_along(samples)) {
  mu_s <- if (cond[s]=="K562") mu_gene*fold_K562 else mu_gene       # element-wise (no ifelse trap)
  mu_tx <- rep(mu_s, ntx) * prop
  reads <- rnbinom(length(tx), mu = pmax(mu_tx,1e-6), size = 15)
  efflen <- pmax(txlen-150, 50); tpm <- reads/efflen; tpm <- if (sum(tpm)>0) tpm/sum(tpm)*1e6 else tpm
  qd <- file.path("results/arm1/salmon", samples[s]); dir.create(qd, recursive=TRUE, showWarnings=FALSE)
  write.table(data.frame(Name=tx, Length=txlen, EffectiveLength=efflen, TPM=round(tpm,4), NumReads=reads),
              file.path(qd,"quant.sf"), sep="\t", quote=FALSE, row.names=FALSE)
}
writeLines(c("conditions:", "  group_A: A549", "  group_B: K562",
             "reference:", "  tx2gene: data/reference/tx2gene.tsv", "threads: 2"), "config/config.yaml")
st <- data.frame(sample_id=samples, condition=cond, platform="illumina", library="short")
write.table(st, "config/samples.tsv", sep="\t", quote=FALSE, row.names=FALSE)
gt <- data.table(gene=gid, is_de=abs(lfc_true)>1, lfc_true=lfc_true)

# ---- run the SHIPPED analysis script ---------------------------------------
log <- system2("Rscript", DE, stdout=TRUE, stderr=TRUE)
if (!file.exists("results/arm1/de/merged_results.tsv")) { cat(log, sep="\n"); stop("deseq2_edger.R produced no output") }

# ---- score -----------------------------------------------------------------
res <- fread("results/arm1/de/merged_results.tsv")
res[, called := !is.na(padj_deseq) & padj_deseq<0.05 & abs(log2FC_deseq)>1]
d <- merge(res, gt, by="gene", all.x=TRUE)
TP <- d[called & is_de, .N]; FP <- d[called & !is_de, .N]; FN <- d[!called & is_de, .N]
prec <- TP/(TP+FP); rec <- TP/(TP+FN)
sign_ok <- d[called & is_de, mean(sign(log2FC_deseq)==sign(lfc_true))]
pass <- prec>=0.90 && rec>=0.80 && sign_ok>=0.99
cat(sprintf("ARM1 DE  precision=%.2f recall=%.2f sign=%.2f  ->  %s\n",
            prec, rec, sign_ok, if(pass)"PASS" else "FAIL"))
writeLines(sprintf('{"arm":"arm1_de","precision":%.3f,"recall":%.3f,"sign":%.3f,"pass":%s}',
                   prec, rec, sign_ok, tolower(as.character(pass))), "metrics.json")
quit(status = if (pass) 0 else 1)
