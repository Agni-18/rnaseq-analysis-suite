# ---- Quantification: Salmon (primary) + STAR (demonstration) ---------------
rule salmon_index:
    input:  txome=config["reference"]["txome_fa"]
    output: directory("results/arm1/salmon_index")
    threads: 4
    shell:  "salmon index -t {input.txome} -i {output} -k 31 -p {threads}"

rule salmon_quant:
    input:  idx="results/arm1/salmon_index",
            fq="results/arm1/trimmed/{sample}.fastq.gz"
    output: "results/arm1/salmon/{sample}/quant.sf"
    threads: 4
    shell:
        "salmon quant -i {input.idx} -l A -r {input.fq} "
        "--validateMappings --gcBias -p {threads} "
        "-o results/arm1/salmon/{wildcards.sample}"

rule star_index:
    input:  fa="data/reference/genome.subset.fa",
            gtf="data/reference/annotation.subset.gtf"
    output: directory("results/arm1/star_index")
    threads: 4
    shell:
        "STAR --runMode genomeGenerate --genomeDir {output} "
        "--genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} "
        "--genomeSAindexNbases 11 --runThreadN {threads}"
