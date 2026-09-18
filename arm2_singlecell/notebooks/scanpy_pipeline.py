#!/usr/bin/env python
"""
Arm 2 — scRNA-seq pipeline: QC -> clustering -> marker-based annotation (Scanpy).

Importable functions (used by the notebook and the validation script) plus a CLI.
Steps mirror the scverse best-practices workflow:
  load -> qc -> (doublets) -> normalize/log/HVG/scale/PCA -> neighbors/UMAP/Leiden
  -> marker-score annotation (each cluster labelled by its top-scoring program).
"""
from __future__ import annotations
import argparse, json
import numpy as np
import scanpy as sc


def load(path: str):
    if path.endswith(".h5ad"):
        adata = sc.read_h5ad(path)
    elif path.endswith(".h5"):
        adata = sc.read_10x_h5(path)
    else:
        adata = sc.read(path)
    adata.var_names_make_unique()
    return adata


def qc(adata, min_genes=200, min_cells=3, max_mt_pct=15.0):
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata = adata[adata.obs["pct_counts_mt"] < max_mt_pct].copy()
    return adata


def detect_doublets(adata):
    try:
        sc.pp.scrublet(adata)
        n = int(adata.obs["predicted_doublet"].sum())
        adata = adata[~adata.obs["predicted_doublet"]].copy()
        print(f">> scrublet removed {n} predicted doublets")
    except Exception as e:
        print(f">> scrublet unavailable, skipping doublet removal ({type(e).__name__})")
    return adata


def preprocess(adata, n_hvg=2000, n_pcs=30):
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)
    adata.raw = adata
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs)
    return adata


def cluster(adata, resolution=1.0, n_pcs=30):
    sc.pp.neighbors(adata, n_pcs=n_pcs)
    sc.tl.umap(adata)
    try:
        sc.tl.leiden(adata, resolution=resolution, flavor="igraph",
                     n_iterations=2, directed=False)
    except TypeError:
        sc.tl.leiden(adata, resolution=resolution)
    return adata


def annotate_markers(adata, markers: dict, key="cell_type"):
    ref_vars = adata.raw.var_names if adata.raw is not None else adata.var_names
    names = list(markers)
    for ct in names:
        genes = [g for g in markers[ct] if g in ref_vars]
        sc.tl.score_genes(adata, genes, score_name=f"score_{ct}")
    scores = adata.obs[[f"score_{ct}" for ct in names]].copy()
    scores.columns = names
    cluster_to_type = scores.groupby(adata.obs["leiden"], observed=True).mean().idxmax(axis=1)
    adata.obs[key] = adata.obs["leiden"].map(cluster_to_type).astype("category")
    return adata


def run(input, output=None, resolution=1.0, n_hvg=2000, n_pcs=30,
        markers=None, do_doublets=True, figdir=None):
    adata = load(input)
    adata = qc(adata)
    if do_doublets:
        adata = detect_doublets(adata)
    adata = preprocess(adata, n_hvg=n_hvg, n_pcs=n_pcs)
    adata = cluster(adata, resolution=resolution, n_pcs=n_pcs)
    if markers:
        adata = annotate_markers(adata, markers)
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
    if figdir:
        import os
        os.makedirs(figdir, exist_ok=True)
        sc.settings.figdir = figdir
        color = ["leiden"] + (["cell_type"] if "cell_type" in adata.obs else [])
        sc.pl.umap(adata, color=color, save="_clusters.png", show=False)
    if output:
        adata.write(output)
    return adata


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Arm 2 scRNA-seq pipeline")
    ap.add_argument("--input", required=True, help="10x .h5 or .h5ad")
    ap.add_argument("--output", default=None, help="output .h5ad")
    ap.add_argument("--markers", default=None, help="JSON {cell_type: [genes...]}")
    ap.add_argument("--resolution", type=float, default=1.0)
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--n-pcs", type=int, default=30)
    ap.add_argument("--no-doublets", action="store_true")
    ap.add_argument("--figdir", default=None)
    a = ap.parse_args()
    mk = json.load(open(a.markers)) if a.markers else None
    run(a.input, output=a.output, resolution=a.resolution, n_hvg=a.n_hvg,
        n_pcs=a.n_pcs, markers=mk, do_doublets=not a.no_doublets, figdir=a.figdir)
