#!/usr/bin/env python
"""Validate the Arm 2 pipeline on a SIMULATED scRNA-seq dataset with known labels.

Realistic hierarchical structure: 2 lineages x 3 sub-types sharing part of their
marker program (like T-cell subsets), a deliberately hard rare subtype (T_reg),
mito genes for QC, and Poisson-Gamma (negative-binomial) counts. Reports
Leiden-vs-truth agreement (ARI/NMI) and marker-annotation accuracy, and writes a
3-panel UMAP (truth / Leiden / predicted) to results/arm2/validation/.

Usage:  python arm2_singlecell/notebooks/validate_arm2.py
"""
import os, sys, json
import numpy as np, scanpy as sc, anndata as ad
from scipy import sparse
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scanpy_pipeline as sp

OUTDIR = "results/arm2/validation"
os.makedirs(OUTDIR, exist_ok=True)
rng = np.random.default_rng(1)

N, G_norm, G_mt = 3500, 1980, 20
lineages = {"L0": ["T_CD4", "T_CD8", "T_reg"], "L1": ["Mono", "DC", "NK"]}
types = [t for ts in lineages.values() for t in ts]; K = len(types)
gene_names = [f"g{j}" for j in range(G_norm)] + [f"MT-{j}" for j in range(G_mt)]
G = G_norm + G_mt
lin_of = {t: lin for lin, ts in lineages.items() for t in ts}

probs = np.array([0.24, 0.20, 0.06, 0.28, 0.04, 0.18]); probs = probs / probs.sum()
labels = rng.choice(K, size=N, p=probs)

base = np.exp(rng.normal(0.2, 1.0, G)); base[G_norm:] *= 0.5
avail = list(range(G_norm)); rng.shuffle(avail)
lineage_prog = {lin: [avail.pop() for _ in range(15)] for lin in lineages}
type_prog = {t: [avail.pop() for _ in range(25)] for t in types}

type_mean = np.tile(base, (K, 1))
for ki, t in enumerate(types):
    type_mean[ki, lineage_prog[lin_of[t]]] *= rng.uniform(1.7, 2.1, size=15)
    if t == "T_reg":                                  # deliberately hard subtype
        type_mean[ki, type_prog[t][:8]] *= rng.uniform(2.0, 2.6, size=8)
    else:
        type_mean[ki, type_prog[t]]     *= rng.uniform(3.0, 6.0, size=25)

X = np.zeros((N, G), dtype=np.int64)
libsize = rng.lognormal(np.log(3000), 0.30, N); shape = 5.0
for i in range(N):
    mu = type_mean[labels[i]].copy(); mu = mu / mu.sum() * libsize[i]
    X[i] = rng.poisson(rng.gamma(shape, mu / shape))

adata = ad.AnnData(X=sparse.csr_matrix(X))
adata.var_names = gene_names
adata.obs_names = [f"cell{i}" for i in range(N)]
adata.obs["true_type"] = [types[l] for l in labels]
markers = {t: [gene_names[j] for j in (lineage_prog[lin_of[t]] + type_prog[t])] for t in types}
print(f">> simulated {N} cells x {G} genes; {K} types in 2 lineages")

adata = sp.qc(adata)
adata = sp.preprocess(adata, n_hvg=2000, n_pcs=30)
adata = sp.cluster(adata, resolution=1.0, n_pcs=30)
adata = sp.annotate_markers(adata, markers)

true = adata.obs["true_type"].values; leiden = adata.obs["leiden"].values; pred = adata.obs["cell_type"].values
ari = adjusted_rand_score(true, leiden); nmi = normalized_mutual_info_score(true, leiden)
acc = float((pred == true).mean())
print(f">> cells after QC: {adata.n_obs} | Leiden clusters: {adata.obs['leiden'].nunique()} (true types: {K})")
print(f">> ARI={ari:.3f}  NMI={nmi:.3f}  annotation_accuracy={acc:.3f}")

sc.settings.figdir = OUTDIR
adata.obs["true_type"] = adata.obs["true_type"].astype("category")
sc.pl.umap(adata, color=["true_type", "leiden", "cell_type"], save="_validation.png", show=False, wspace=0.4)
json.dump({"ari": ari, "nmi": nmi, "accuracy": acc, "n_cells": int(adata.n_obs),
           "n_clusters": int(adata.obs['leiden'].nunique()), "n_types": K},
          open(os.path.join(OUTDIR, "metrics.json"), "w"), indent=2)
print(f">> wrote {OUTDIR}/metrics.json and umap_validation.png")
