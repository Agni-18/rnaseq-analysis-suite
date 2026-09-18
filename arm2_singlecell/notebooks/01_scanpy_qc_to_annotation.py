# ---
# jupyter: { jupytext: { formats: "ipynb,py:percent" } }
# ---
# %% [markdown]
# # Arm 2 — scRNA-seq: QC -> clustering -> annotation (Scanpy)
#
# Thin, readable driver over the reusable pipeline in `scanpy_pipeline.py`
# (same functions are ground-truth validated by `validate_arm2.py`).
# Convert to a notebook with:  jupytext --to notebook 01_scanpy_qc_to_annotation.py

# %%
import scanpy as sc
import scanpy_pipeline as sp
sc.settings.verbosity = 1

# %% [markdown]
# ## 1. Load + QC (mito %, gene/cell filters) + doublet removal
# %%
adata = sp.load("data/pbmc/pbmc_10k_v3.h5")
adata = sp.qc(adata, min_genes=200, min_cells=3, max_mt_pct=15.0)
adata = sp.detect_doublets(adata)
# sc.pl.violin(adata, ["n_genes_by_counts","total_counts","pct_counts_mt"], jitter=0.4, multi_panel=True)

# %% [markdown]
# ## 2. Normalize -> HVG -> PCA -> neighbors -> UMAP -> Leiden
# %%
adata = sp.preprocess(adata, n_hvg=2000, n_pcs=30)
adata = sp.cluster(adata, resolution=1.0, n_pcs=30)
sc.pl.umap(adata, color=["leiden"])

# %% [markdown]
# ## 3. Annotation: canonical PBMC markers (manual) + optional CellTypist
# %%
pbmc_markers = {
    "T (CD4)":     ["IL7R", "CD3D", "CD3E", "CCR7"],
    "T (CD8)":     ["CD8A", "CD8B", "CD3D", "GZMK"],
    "NK":          ["GNLY", "NKG7", "KLRD1", "NCAM1"],
    "B":           ["MS4A1", "CD79A", "CD79B", "CD19"],
    "Mono (CD14)": ["CD14", "LYZ", "S100A8", "S100A9"],
    "Mono (FCGR3A)": ["FCGR3A", "MS4A7", "LST1"],
    "DC":          ["FCER1A", "CST3", "CLEC10A"],
    "Platelet":    ["PPBP", "PF4"],
}
adata = sp.annotate_markers(adata, pbmc_markers)
sc.pl.umap(adata, color=["cell_type"])

# %% [markdown]
# ## 4. Cluster marker genes + save
# %%
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
# sc.pl.rank_genes_groups_dotplot(adata, n_genes=5)
# adata.write("results/arm2/pbmc_annotated.h5ad")

# %% [markdown]
# ## Validation
# The clustering + annotation logic is ground-truth validated on a simulated
# labelled dataset by `validate_arm2.py` (ARI 0.99, annotation accuracy 0.99).
