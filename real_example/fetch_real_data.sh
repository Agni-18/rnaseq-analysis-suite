#!/usr/bin/env bash
# Fetch the real fission-yeast RNA-seq dataset (Data Carpentry / Bahler lab) used
# by the Arm 1 real example. These GitHub-hosted files are reachable without the
# SG-NEx/GENCODE hosts. (Arm 2's real PBMC data ships offline inside Scanpy.)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
DATA="$HERE/data"; mkdir -p "$DATA"
BASE="https://raw.githubusercontent.com/tavareshugo/data-carpentry-rnaseq/master/data"
for f in counts_raw.csv sample_info.csv; do
  echo ">> fetching $f"
  curl -fsSL -o "$DATA/$f" "$BASE/$f"
  [ -s "$DATA/$f" ] || { echo "!! $f is empty"; exit 1; }
done
echo ">> real data in $DATA:  $(wc -l < "$DATA/counts_raw.csv") count rows, $(( $(wc -l < "$DATA/sample_info.csv") - 1 )) samples"
