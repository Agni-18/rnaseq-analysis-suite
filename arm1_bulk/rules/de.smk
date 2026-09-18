# ---- Differential expression + functional (handoff to R) --------------------
rule tx2gene:
    input:  txome=config["reference"]["txome_fa"]
    output: "data/reference/tx2gene.tsv"
    shell:  "bash scripts/make_tx2gene.sh {input.txome} {output}"

rule differential_expression:
    input:  quant=expand("results/arm1/salmon/{s}/quant.sf", s=ILLUMINA),
            t2g="data/reference/tx2gene.tsv"
    output: "results/arm1/de/deseq2_vs_edger_concordance.png",
            rnk="results/arm1/de/ranked_genes.rnk"
    conda:  "../../environments/env_bulk.yml"
    shell:  "Rscript arm1_bulk/R/deseq2_edger.R"

rule functional:
    input:  "results/arm1/de/ranked_genes.rnk"
    output: "results/arm1/functional/gsea_dotplot.png"
    conda:  "../../environments/env_bulk.yml"
    shell:  "Rscript arm1_bulk/R/functional.R"
