#!/usr/bin/env bash
# Fetch GRCh38 reference (genome, GTF, transcriptome) from GENCODE and, for a
# laptop, subset to the chromosomes in config/config.yaml (default chr20-22).
set -euo pipefail
mkdir -p data/reference && cd data/reference
GENCODE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46"
[ -f genome.fa ]        || { wget -q "${GENCODE}/GRCh38.primary_assembly.genome.fa.gz"     -O genome.fa.gz && pigz -d genome.fa.gz; }
[ -f annotation.gtf ]   || { wget -q "${GENCODE}/gencode.v46.primary_assembly.annotation.gtf.gz" -O annotation.gtf.gz && pigz -d annotation.gtf.gz; }
[ -f transcriptome.fa ] || { wget -q "${GENCODE}/gencode.v46.transcripts.fa.gz"             -O transcriptome.fa.gz && pigz -d transcriptome.fa.gz; }
CHROMS="${1:-chr20 chr21 chr22}"
if [ -n "${CHROMS}" ]; then
  echo ">> Subsetting genome + GTF to: ${CHROMS}"
  samtools faidx genome.fa
  samtools faidx genome.fa ${CHROMS} > genome.subset.fa
  grep -E "^($(echo ${CHROMS} | tr ' ' '|'))\b" annotation.gtf > annotation.subset.gtf
fi
echo ">> References ready in data/reference/"
