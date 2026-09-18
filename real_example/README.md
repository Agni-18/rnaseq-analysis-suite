# Real-data example (reproducible without the SG-NEx / GENCODE / 10x hosts)

The main pipeline targets **SG-NEx** (bulk + Nanopore) and **10x PBMC**. Those
hosts — and the conda channels — are not always reachable (e.g. locked-down CI).
This folder runs the **same analysis code** on real public datasets that are
reachable, so the methodology is demonstrable end-to-end anywhere.

```
make real-example      # from the repo root: fetch + run both arms
```
or step by step:
```
bash real_example/fetch_real_data.sh     # real fission-yeast counts (GitHub)
Rscript real_example/real_arm1_de.R      # DESeq2 + edgeR + concordance
python real_example/real_arm2_pbmc.py    # Scanpy clustering + annotation
```
Outputs land in `real_example/results/`. Committed reference figures (what you
should get) are in `real_example/example_outputs/`.

## Arm 1 — bulk DE, real *S. pombe* wt-vs-mutant stress timecourse
Same core as `arm1_bulk/R/deseq2_edger.R` (DESeq2 `lfcShrink` + edgeR TMM/glmQLF
+ concordance), count-matrix entry point (estimated counts rounded to integers).
Two contrasts:
- **Stress response, wt 0 vs 180 min (3v3):** DESeq2 32, edgeR 76, shared 32,
  Jaccard 0.42, r(log2FC) 0.93 at padj<0.05 & |log2FC|>1.
- **Mutant vs wild-type, all timepoints (~minute+strain, 18v18):** DESeq2 125,
  edgeR 81, shared 80, Jaccard 0.63, r 0.87 at padj<0.05. Highly significant but
  small-magnitude — a reminder to match the effect-size cutoff to the contrast.
The PCA shows PC1 (61%) is the time axis; strain is a subtle secondary effect.

## Arm 2 — scRNA-seq, real 10x PBMC
`pbmc68k_reduced` ships inside Scanpy (no download) with the study's cell-type
labels. Runs the shipped `scanpy_pipeline.py` clustering + marker annotation:
**ARI 0.40, NMI 0.61** vs the study labels. Major lineages (B / Monocyte /
Dendritic) resolve cleanly; fine T-cell subsets overlap — the real-world
counterpart to the simulated validation (ARI 0.99, which proves code correctness).

## Note on the fold-change fan
DESeq2 log2FCs are shrunk (`lfcShrink`) while edgeR's glmQLF log2FCs are not, so
edgeR's are slightly larger in magnitude — the fan around the diagonal in the
concordance plots is expected.
