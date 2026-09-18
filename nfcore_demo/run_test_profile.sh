#!/usr/bin/env bash
# Flourish — run nf-core/rnaseq on its bundled test profile (tiny data, laptop-safe).
# Signals fluency with the de-facto industry-standard pipeline these shops run.
set -euo pipefail
nextflow run nf-core/rnaseq -profile test,docker --outdir results/nfcore_rnaseq
echo ">> nf-core/rnaseq test complete. Compare its MultiQC to Arm 1's."
