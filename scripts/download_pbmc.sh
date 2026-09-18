#!/usr/bin/env bash
# Fetch a 10x human PBMC dataset for Arm 2 (single-cell).
# Simplest path is Scanpy's built-in 3k PBMC:
#   python -c "import scanpy as sc; sc.datasets.pbmc3k()"
# For a larger ~10k set, pull a 10x filtered feature-barcode matrix:
set -euo pipefail
mkdir -p data/pbmc && cd data/pbmc
URL="https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5"
[ -f pbmc_10k_v3.h5 ] || wget -q "$URL" -O pbmc_10k_v3.h5
echo ">> PBMC 10k matrix at data/pbmc/pbmc_10k_v3.h5"
