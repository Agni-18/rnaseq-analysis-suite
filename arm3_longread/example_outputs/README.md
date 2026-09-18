# Arm 3 — example outputs (DTU)

Figures produced by `arm3_longread/R/dtu.R` (DRIMSeq + two-stage stageR) on a
**simulated Bambu-style transcript-count matrix** with known isoform-switch
ground truth: 800 genes (2-4 transcripts each), 5 vs 5 samples, 120 spiked-in
DTU genes. DTU genes hold **total** gene expression matched across conditions
and change only isoform **proportions** — so this tests usage, not DGE.

Validation (gene-level DTU vs ground truth): **precision 0.91, recall 0.97**
(deterministic; reproduce with `make validate`).

`top_switch_gene_proportions.png` is the signature long-read result: an isoform
whose usage flips between conditions — a change a gene-level short-read test
(Arm 1) cannot see. `dtu_gene_pvalues.png` is the gene-level p-value histogram.

Note: DRIMSeq's gene-level LRT is mildly anti-conservative at small n, so
stageR's 5% OFDR is a nominal target. On real SG-NEx A549-vs-K562 data these
regenerate from Bambu counts with the same code (p2-longread env).
