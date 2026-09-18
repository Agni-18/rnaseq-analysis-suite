#!/usr/bin/env python
"""REAL-data example — Arm 2 scRNA-seq on real 10x PBMC (Scanpy's bundled,
offline pbmc68k_reduced, which ships the original study's cell-type labels).
Runs the shipped arm2_singlecell/notebooks/scanpy_pipeline.py clustering +
marker annotation and scores against the study labels.
Usage:  python real_example/real_arm2_pbmc.py
"""
import os, sys, json, warnings; warnings.filterwarnings("ignore")
import scanpy as sc
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "arm2_singlecell", "notebooks"))
import scanpy_pipeline as sp

OUT = os.path.join(HERE, "results"); os.makedirs(OUT, exist_ok=True)
sc.settings.verbosity = 0
adata = sc.datasets.pbmc68k_reduced()          # real 10x PBMC, offline, with 'bulk_labels'
print(f">> real PBMC: {adata.n_obs} cells x {adata.n_vars} genes")

sc.pp.neighbors(adata, n_pcs=30, use_rep="X_pca")
try:
    sc.tl.leiden(adata, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)
except TypeError:
    sc.tl.leiden(adata, resolution=1.0)
if "X_umap" not in adata.obsm:
    sc.tl.umap(adata)

pbmc_markers = {
    "T cell":     ["CD3D","CD3E","IL7R","CD2"],
    "CD8 T / NK": ["CD8A","GZMB","NKG7","GNLY","KLRD1"],
    "B cell":     ["CD79A","CD79B","MS4A1","CD19"],
    "Monocyte":   ["LYZ","S100A8","S100A9","CST3","FCN1"],
    "Dendritic":  ["FCER1A","CLEC10A"],
}
avail = set(adata.raw.var_names if adata.raw is not None else adata.var_names)
pbmc_markers = {k: [g for g in v if g in avail] for k, v in pbmc_markers.items()}
adata = sp.annotate_markers(adata, pbmc_markers)

truth, leiden = adata.obs["bulk_labels"].values, adata.obs["leiden"].values
ari = adjusted_rand_score(truth, leiden); nmi = normalized_mutual_info_score(truth, leiden)
print(f">> REAL DATA  ARI(Leiden, study labels)={ari:.3f}  NMI={nmi:.3f}  clusters={adata.obs['leiden'].nunique()}")

sc.settings.figdir = OUT
sc.pl.umap(adata, color=["bulk_labels","leiden","cell_type"], save="_real_pbmc.png", show=False, wspace=0.5)
json.dump({"dataset":"pbmc68k_reduced (real 10x PBMC)","n_cells":int(adata.n_obs),
           "ari":round(float(ari),3),"nmi":round(float(nmi),3),
           "n_clusters":int(adata.obs['leiden'].nunique())},
          open(os.path.join(OUT,"real_arm2_metrics.json"),"w"), indent=2)
print(">> Arm 2 real example done. UMAP + metrics in real_example/results/")
