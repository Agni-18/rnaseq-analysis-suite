# ---- QC + trimming ----------------------------------------------------------
rule fastqc:
    input:  "data/illumina/{sample}/{sample}.fastq.gz"
    output: "results/arm1/qc/{sample}_fastqc.zip"
    threads: 2
    shell:  "fastqc -t {threads} -o results/arm1/qc {input}"

rule fastp:
    input:  "data/illumina/{sample}/{sample}.fastq.gz"
    output: trimmed="results/arm1/trimmed/{sample}.fastq.gz",
            html="results/arm1/qc/{sample}_fastp.html"
    threads: 4
    shell:  "fastp -i {input} -o {output.trimmed} -h {output.html} -w {threads}"

rule multiqc:
    input:  expand("results/arm1/qc/{s}_fastqc.zip", s=ILLUMINA)
    output: "results/arm1/multiqc_report.html"
    shell:  "multiqc results/arm1/qc -n multiqc_report.html -o results/arm1"
