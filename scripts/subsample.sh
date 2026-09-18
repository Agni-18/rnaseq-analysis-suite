#!/usr/bin/env bash
# Downsample FASTQs with seqtk so pipelines run on a laptop.
# Usage: scripts/subsample.sh <in.fastq.gz> <out.fastq.gz> <n_reads> [seed]
set -euo pipefail
IN="$1"; OUT="$2"; N="${3:-500000}"; SEED="${4:-42}"
[ "$N" -eq 0 ] && { cp "$IN" "$OUT"; echo ">> downsampling disabled, copied $IN"; exit 0; }
seqtk sample -s"${SEED}" "$IN" "$N" | pigz > "$OUT"
echo ">> ${IN} -> ${OUT} (${N} reads, seed ${SEED})"
