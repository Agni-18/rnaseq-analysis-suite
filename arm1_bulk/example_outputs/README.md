# Arm 1 — example outputs

Figures produced by `arm1_bulk/R/deseq2_edger.R` on a **simulated Salmon dataset**
(6 samples, 3 vs 3; 1500 genes; 200 spiked-in DE genes with known ground truth).

Validation (DESeq2 vs ground truth): **precision 0.98, recall 0.92, sign 100%**.
Engine concordance: **Jaccard 0.97, r(log2FC) 0.998**.

On real SG-NEx A549-vs-K562 data these plots regenerate with the same code.
(The validation sandbox lacked `apeglm`, so shrinkage fell back to the normal
prior; the production `environments/env_bulk.yml` pins `apeglm`, which the
concordance-plot axis label assumes.)
