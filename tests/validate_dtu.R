#!/usr/bin/env Rscript
# =============================================================================
# Arm 3 validation: simulate transcript counts with known isoform-switch ground
# truth, run the SHIPPED arm3_longread/R/dtu.R, and score gene-level DTU
# precision/recall. Exits non-zero if thresholds aren't met.
# Usage:  Rscript tests/validate_dtu.R
# =============================================================================
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = FALSE)
self <- sub("^--file=", "", args[grep("^--file=", args)])
repo <- normalizePath(file.path(dirname(self), ".."))
DTU  <- file.path(repo, "arm3_longread", "R", "dtu.R")
work <- file.path(repo, "tests", ".work", "arm3")
unlink(work, recursive = TRUE); dir.create(file.path(work, "config"), recursive = TRUE)
setwd(work)

# ---- simulate (matched to the documented Arm 3 validation) ------------------
set.seed(7)
G <- 800L; N_DTU <- 120L; conc <- 200; reps <- 1:5
samples <- c(paste0("SGNex_A549_directcDNA_replicate", reps, "_run1"),
             paste0("SGNex_K562_directcDNA_replicate", reps, "_run1"))
cond <- c(rep("A549",5), rep("K562",5))
rdirichlet <- function(a){ x <- rgamma(length(a), a, 1); x/sum(x) }
ntx <- sample(2:4, G, TRUE); gene_id <- sprintf("G%04d", 1:G)
gene_mu <- exp(rnorm(G, 6.8, 0.7))
dtu_idx <- sample(G, N_DTU); is_dtu <- logical(G); is_dtu[dtu_idx] <- TRUE
shift <- numeric(G); shift[dtu_idx] <- runif(N_DTU, 0.6, 1.0)
feat <- unlist(lapply(seq_len(G), function(i) paste0(gene_id[i], "_t", seq_len(ntx[i]))))
feat_gene <- rep(gene_id, ntx)
propA <- lapply(seq_len(G), function(i){ r <- sort(runif(ntx[i], .4, 3))^2; r/sum(r) })
propB <- lapply(seq_len(G), function(i){ if(!is_dtu[i]) return(propA[[i]]); s<-shift[i]; pr<-(1-s)*propA[[i]]+s*rev(propA[[i]]); pr/sum(pr) })
cmat <- matrix(0L, nrow=length(feat), ncol=length(samples), dimnames=list(feat, samples))
row0 <- 1L
for(i in seq_len(G)){
  rows <- row0:(row0+ntx[i]-1L); row0 <- row0+ntx[i]
  for(s in seq_along(samples)){
    base <- if(cond[s]=="K562") propB[[i]] else propA[[i]]
    tot <- rnbinom(1, mu=gene_mu[i], size=20)
    if(tot>0) cmat[rows,s] <- as.integer(rmultinom(1, tot, rdirichlet(base*conc)))
  }
}
dir.create("results/arm3/counts", recursive=TRUE, showWarnings=FALSE)
fwrite(data.table(feature_id=feat, gene_id=feat_gene, as.data.table(cmat)),
       "results/arm3/counts/transcript_counts.tsv", sep="\t")
writeLines(c("conditions:", "  group_A: A549", "  group_B: K562", "threads: 2"), "config/config.yaml")
st <- data.frame(sample_id=samples, condition=cond, platform="nanopore", library="cdna")
write.table(st, "config/samples.tsv", sep="\t", quote=FALSE, row.names=FALSE)
gt <- data.table(gene_id=gene_id, is_dtu=is_dtu)

# ---- run the SHIPPED analysis script ---------------------------------------
log <- system2("Rscript", DTU, stdout=TRUE, stderr=TRUE)
if (!file.exists("results/arm3/dtu/dtu_results.tsv")) { cat(log, sep="\n"); stop("dtu.R produced no output") }

# ---- score -----------------------------------------------------------------
res <- fread("results/arm3/dtu/dtu_results.tsv")
called <- unique(res[padj_gene_screen < 0.05, gene_id]); gt[, called := gene_id %in% called]
TP <- gt[is_dtu & called, .N]; FP <- gt[!is_dtu & called, .N]; FN <- gt[is_dtu & !called, .N]
prec <- TP/(TP+FP); rec <- TP/(TP+FN)
pass <- prec>=0.80 && rec>=0.90
cat(sprintf("ARM3 DTU precision=%.2f recall=%.2f  ->  %s\n", prec, rec, if(pass)"PASS" else "FAIL"))
writeLines(sprintf('{"arm":"arm3_dtu","precision":%.3f,"recall":%.3f,"pass":%s}',
                   prec, rec, tolower(as.character(pass))), "metrics.json")
quit(status = if (pass) 0 else 1)
