# Validation

Each arm's analysis code is exercised on a simulated dataset with known
ground truth by a script under `tests/` (Arm 2's lives with its module).
Regenerate this file any time with `make validate`.

| Arm | Result vs ground truth | Thresholds | Status |
|-----|------------------------|------------|--------|
| Arm 1 — bulk DE (DESeq2+edgeR) | precision 0.975, recall 0.924, sign 1.0 | precision>=0.90, recall>=0.80, sign>=0.99 | PASS |
| Arm 2 — scRNA-seq (Leiden+markers) | ARI 0.992, NMI 0.983, accuracy 0.995 | ARI>=0.90, accuracy>=0.90 | PASS |
| Arm 3 — long-read DTU (DRIMSeq+stageR) | precision 0.914, recall 0.975 | precision>=0.80, recall>=0.90 | PASS |

**Overall: PASS**
