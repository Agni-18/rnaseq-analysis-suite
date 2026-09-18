#!/usr/bin/env Rscript
# Arm 2 (parallel) — compact Seurat v5 pass mirroring the Scanpy pipeline, to
# show cross-tool fluency and let the two label sets be compared.
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })

so <- Read10X_h5("data/pbmc/pbmc_10k_v3.h5") |>
  CreateSeuratObject(min.cells = 3, min.features = 200)
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so <- subset(so, subset = nFeature_RNA > 200 & percent.mt < 15)

so <- SCTransform(so, verbose = FALSE) |>
  RunPCA(verbose = FALSE) |>
  FindNeighbors(dims = 1:30) |>
  FindClusters(resolution = 1.0) |>
  RunUMAP(dims = 1:30)

markers <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# SingleR or manual annotation here, then cross-tabulate Seurat clusters vs the
# Scanpy Leiden labels (export both to compare concordance).
dir.create("results/arm2", recursive = TRUE, showWarnings = FALSE)
saveRDS(so, "results/arm2/seurat_object.rds")
write.csv(markers, "results/arm2/seurat_cluster_markers.csv", row.names = FALSE)
message(">> Seurat parallel done: UMAP + cluster markers written to results/arm2/")
