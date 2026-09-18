#!/usr/bin/env bash
# Download SG-NEx data (matched Illumina + Nanopore, human cell lines) for Arms 1 & 3.
# Public bucket, no AWS account needed (--no-sign-request).
# Data source: registry.opendata.aws/sg-nex-data  (GRCh38-aligned, GoekeLab)
# LIST first, then edit the sample lists below with the aliases you want.
set -euo pipefail
BUCKET="s3://sg-nex-data/data"
mkdir -p data/illumina data/nanopore data/reference
echo ">> Available Illumina samples (short-read):"
aws s3 ls --no-sign-request "${BUCKET}/sequencing_data_illumina/fastq/" | grep -Ei 'A549|K562' || true
echo ">> Available Nanopore samples (long-read cDNA):"
aws s3 ls --no-sign-request "${BUCKET}/sequencing_data_ont/fastq/" | grep -Ei 'A549|K562' | grep -i cdna || true
ILLUMINA_SAMPLES=(
  # SGNex_A549_Illumina_replicate1_run1
  # SGNex_K562_Illumina_replicate1_run1
)
NANOPORE_SAMPLES=(
  # SGNex_A549_directcDNA_replicate1_run1
  # SGNex_K562_directcDNA_replicate1_run1
)
for s in "${ILLUMINA_SAMPLES[@]}"; do
  echo ">> Illumina: $s"
  aws s3 sync --no-sign-request "${BUCKET}/sequencing_data_illumina/fastq/${s}" "data/illumina/${s}"
done
for s in "${NANOPORE_SAMPLES[@]}"; do
  echo ">> Nanopore: $s"
  aws s3 sync --no-sign-request "${BUCKET}/sequencing_data_ont/fastq/${s}" "data/nanopore/${s}"
done
echo ">> Done. Cite: SG-NEx accessed on $(date +%F) at registry.opendata.aws/sg-nex-data"
