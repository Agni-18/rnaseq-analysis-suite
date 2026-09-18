# ---- Isoform discovery/quant -> DTU -----------------------------------------
# Bambu (R) is the primary discovery/quant engine; it writes the transcript-count
# contract (feature_id, gene_id, <samples>) that dtu.R consumes.
rule bambu:
    input:  bams=expand("results/arm3/bam/{s}.bam", s=NANOPORE),
            ref="data/reference/genome.subset.fa",
            gtf="data/reference/annotation.subset.gtf"
    output: "results/arm3/counts/transcript_counts.tsv"
    conda:  "../../environments/env_longread.yml"
    shell:  "Rscript arm3_longread/R/run_bambu.R"

# IsoQuant is an alternative discovery/quant engine (run alongside for comparison).
rule isoquant:
    input:  bams=expand("results/arm3/bam/{s}.bam", s=NANOPORE),
            ref="data/reference/genome.subset.fa",
            gtf="data/reference/annotation.subset.gtf"
    output: "results/arm3/isoquant/OUT/OUT.transcript_models.gtf"
    threads: 4
    shell:
        "isoquant.py --reference {input.ref} --genedb {input.gtf} "
        "--bam {input.bams} --data_type nanopore -t {threads} "
        "-o results/arm3/isoquant"

# Differential transcript usage: DRIMSeq + two-stage stageR.
rule dtu:
    input:  "results/arm3/counts/transcript_counts.tsv"
    output: "results/arm3/dtu/dtu_results.tsv",
            "results/arm3/dtu/top_switch_gene_proportions.png"
    conda:  "../../environments/env_longread.yml"
    shell:  "Rscript arm3_longread/R/dtu.R"
