# ---- QC + splice-aware alignment (minimap2) ---------------------------------
rule nanoplot:
    input:  "data/nanopore/{sample}/{sample}.fastq.gz"
    output: "results/arm3/qc/{sample}_NanoPlot/NanoStats.txt"
    threads: 2
    shell:  "NanoPlot -t {threads} --fastq {input} -o results/arm3/qc/{wildcards.sample}_NanoPlot"

rule minimap2_splice:
    input:  fq="data/nanopore/{sample}/{sample}.fastq.gz",
            ref="data/reference/genome.subset.fa"
    output: "results/arm3/bam/{sample}.bam"
    threads: 4
    # -ax splice for spliced cDNA long reads; add -uf -k14 for direct-RNA.
    shell:
        "minimap2 -ax splice -t {threads} {input.ref} {input.fq} "
        "| samtools sort -@ {threads} -o {output} && samtools index {output}"
