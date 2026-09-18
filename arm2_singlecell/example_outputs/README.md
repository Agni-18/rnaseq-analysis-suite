# Arm 2 — example outputs (single-cell)

`umap_validation.png` is produced by `arm2_singlecell/notebooks/validate_arm2.py`,
which runs the exact pipeline functions from `scanpy_pipeline.py` on a **simulated
labelled dataset**: 3500 cells, 6 cell types in 2 lineages that share part of
their marker program (like T-cell subsets), a deliberately hard rare subtype
(T_reg), mito genes for QC, and Poisson-Gamma counts.

Validation vs known labels: **ARI 0.99, NMI 0.98, annotation accuracy 0.99**
(see `metrics.json`). The three UMAP panels are truth / Leiden / marker-predicted.

Caveat: transcriptionally distinct populations separate near-perfectly under
standard Leiden, so this favourable number reflects distinct-ish simulated types.
Real data (ambient RNA, batch effects, continuous states, fine subsets) is harder
and is where resolution tuning and integration matter. The point here is that the
QC -> cluster -> annotate logic recovers known structure correctly.
