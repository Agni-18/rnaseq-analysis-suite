#!/usr/bin/env bash
# Build tx2gene.tsv (transcript_id  gene_id  gene_symbol) from the GENCODE
# transcriptome FASTA used to build the Salmon index. Pipe-delimited headers:
# >ENST..|ENSG..|OTTHUMG..|OTTHUMT..|txname|GENE_SYMBOL|len|...
set -euo pipefail
FA="${1:-data/reference/transcriptome.fa}"
OUT="${2:-data/reference/tx2gene.tsv}"
grep '^>' "$FA" | sed 's/^>//' | awk -F'|' 'BEGIN{OFS="\t"}{print $1,$2,$6}' > "$OUT"
echo ">> wrote $(wc -l < "$OUT") rows to $OUT"
