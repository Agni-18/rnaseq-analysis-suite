# Project 2 — RNA-seq / Functional Genomics + Single-Cell (+ Long-Read Isoform)

**Status:** Draft for sign-off
**Author:** (you)
**Target employers:** Genotypic Technology, MedGenome, Strand Life Sciences, Eurofins

---

## 1. One-line purpose

A three-arm transcriptomics portfolio project that demonstrates end-to-end range: a **bulk RNA-seq differential-expression pipeline in R**, a **single-cell RNA-seq analysis in Python**, and a **Nanopore long-read isoform component** as a technical flourish — packaged with a real workflow manager to signal production readiness.

---

## 2. Why this project, and who it speaks to

| Arm | Skills shown | Primary target fit |
|-----|--------------|--------------------|
| Bulk RNA-seq (R) | STAR/Salmon, DESeq2/edgeR, GO/KEGG/GSEA, statistics | Genotypic, Eurofins, Strand |
| scRNA-seq (Python) | Scanpy/Seurat, QC, clustering, annotation, markers | MedGenome, Strand |
| Long-read isoform (flourish) | minimap2, IsoQuant/Bambu, DTU, novel-isoform discovery | Genotypic (ONT service core) |
| Packaging | Snakemake, nf-core fluency, conda, reproducibility | All four (service/production shops) |

The **bulk-in-R / single-cell-in-Python** split deliberately shows two toolchains. The **shared biology between the bulk and long-read arms** (see §4) lets us demonstrate short-vs-long-read complementarity — the exact value proposition a long-read service provider like Genotypic sells.

---

## 3. Scope & constraints

- **Compute target:** a single laptop with ≤16 GB RAM. This is treated as a *design driver*, not an afterthought.
- **Strategy:** lightweight-first tools (Salmon, minimap2), a **chromosome subset (chr20–22)** for the memory-heavy demonstrations (STAR alignment), and **downsampled reads** everywhere. Every pipeline is written to be depth-agnostic and genome-scalable, documented as "scales to full genome on HPC."
- **Honesty framing:** all subsetting/downsampling is disclosed in the README as a compute-scoping choice, not hidden. This reads as good practice to a reviewer.

---

## 4. Data sources

### Arms 1 & 3 — shared biology: SG-NEx
The **Singapore Nanopore Expression (SG-NEx)** resource provides **matched Illumina short-read and Oxford Nanopore long-read** data for the *same* human cell lines (e.g. A549, K562, MCF7, HepG2), plus synthetic spike-ins (sequins/SIRVs) useful for sanity checks. Public on AWS Open Data.

- **Bulk arm** uses the Illumina short-read libraries.
- **Long-read arm** uses the ONT libraries (direct cDNA / PCR-cDNA; direct RNA optional).
- Because they are the same samples, we can put a **short-vs-long complementarity figure** in the final report.

Pick **two contrasting cell lines** (e.g. A549 vs K562) so the bulk arm has a real two-group differential-expression contrast with replicates.

### Arm 2 — single-cell: 10x PBMC
A ~5–10k-cell **10x Genomics human PBMC** dataset. Laptop-tractable, canonical immune markers, and **CellTypist** ships ready human-immune reference models for the automated-annotation showcase. Swappable for a small tumor/tissue dataset if we want a stronger disease narrative later.

---

## 5. Architecture overview

```
                    ┌─────────────────────────────┐
                    │        SG-NEx (shared)       │
                    │  Illumina  +  Nanopore ONT   │
                    └───────────┬─────────┬────────┘
                                │         │
              ┌─────────────────▼──┐   ┌──▼──────────────────┐
   ARM 1      │  BULK RNA-seq (R)  │   │  LONG-READ (ARM 3)  │   ARM 3
              │  Snakemake DAG     │   │  Snakemake DAG      │
              └─────────┬──────────┘   └──────────┬──────────┘
                        │                         │
                        └──────────┬──────────────┘
                                   ▼
                      short-vs-long complementarity report

   ARM 2   10x PBMC ──►  scRNA-seq (Python, Scanpy/Seurat notebooks)

   FLOURISH  nf-core/rnaseq  -profile test   (industry-standard signal)
```

Snakemake wraps the two deterministic arms (bulk, long-read). Single-cell stays in reviewable notebooks because that is how exploratory sc analysis is genuinely delivered.

---

## 6. Arm 1 — Bulk RNA-seq (R-centric)

**Goal:** a two-condition differential-expression + functional-enrichment pipeline, with two DE engines cross-checked.

**Steps**
1. **QC:** FastQC → aggregated with MultiQC.
2. **Trimming:** fastp (adapter + quality).
3. **Quantification (primary):** Salmon selective alignment against a transcriptome index (GC + sequence-bias correction), a few GB RAM.
4. **Quantification (demonstration):** STAR on a **chr20–22 subset index** → featureCounts (subread). Shows alignment/BAM/IGV competency without the full-genome memory wall.
5. **Import to R:** tximport (Salmon → gene-level counts).
6. **Differential expression, two engines:**
   - **DESeq2** — median-of-ratios normalization, negative-binomial GLM, Wald test, `apeglm` LFC shrinkage.
   - **edgeR** — TMM normalization, quasi-likelihood F-test.
   - **Concordance analysis** of the two DEG sets (overlap + logFC agreement) as a deliberate statistics showcase.
7. **Functional layer:**
   - GO + KEGG over-representation with **clusterProfiler**.
   - **GSEA** on ranked log2FC with **fgsea** / clusterProfiler.

**Deliverables:** MultiQC report; PCA; sample-correlation heatmap; MA plot; volcano; DESeq2-vs-edgeR concordance figure; GO/KEGG dotplots; GSEA enrichment plots; a results table.

---

## 7. Arm 2 — scRNA-seq (Python-centric)

**Goal:** a standard-of-practice single-cell workflow from raw counts to annotated cell types, primary in Scanpy with a compact Seurat parallel.

**Steps (Scanpy primary)**
1. **QC:** filter on min genes/cell, min cells/gene, mitochondrial %; QC violin plots.
2. **Doublet removal:** Scrublet (or scDblFinder in the Seurat pass).
3. **(Optional) ambient RNA:** SoupX / CellBender note.
4. **Normalize:** `normalize_total` + `log1p`.
5. **Feature selection:** highly variable genes.
6. **Dimensionality reduction:** PCA → neighbors → UMAP.
7. **Clustering:** **Leiden**.
8. **Annotation, two ways:**
   - Manual: `rank_genes_groups` (Wilcoxon) + canonical immune markers.
   - Automated: **CellTypist** (human immune models).
   - Compare manual vs automated as a showcase.
9. **Integration (if multi-sample):** Harmony or scVI (scvi-tools).

**Seurat parallel:** a shorter Seurat v5 pass (LogNormalize/SCTransform → clustering → markers) to prove cross-tool fluency.

**Deliverables:** QC violins; UMAPs colored by cluster / sample / annotation; marker dotplot + heatmap; manual-vs-CellTypist annotation comparison; top-markers table.

---

## 8. Arm 3 — Long-read isoform (the flourish)

**Goal:** show the thing short reads cannot do — resolve full-length isoforms, discover novel ones, and test differential transcript usage.

**Steps** (basecalling assumed pre-done with Dorado; start from FASTQ)
1. **QC:** NanoPlot / pycoQC (read-length, quality, yield).
2. **Subsample:** seqtk to a laptop-tractable depth.
3. **Alignment:** **minimap2** splice-aware (`-ax splice`; direct-RNA settings if used), against the chr20–22 subset.
4. **Isoform discovery + quantification:** **IsoQuant** and/or **Bambu** (reference-guided; current-generation choices for well-annotated human genome).
5. **Differential transcript usage (DTU):** DRIMSeq + DEXSeq + stageR (or limma diffSplice / satuRn).
6. **Complementarity:** pick a gene where long reads resolve isoforms the short-read arm collapses; show both.

**Deliverables:** NanoPlot summary; alignment stats; known-vs-novel isoform counts; an IGV isoform-structure snapshot; DTU / isoform-switch example; the short-vs-long complementarity figure.

**Design note:** reference-guided tools are used deliberately — current benchmarks (LRGASP) find reference-based methods perform best in well-annotated genomes like human, and that longer/more accurate reads matter more than raw depth for transcript models.

---

## 9. Packaging & reproducibility

- **Snakemake** DAG for Arms 1 and 3 (rules, wildcards, config-driven samples). Writing our *own* DAG proves pipeline-*building*, not just pipeline-*running*.
- **Notebooks** (Jupyter) for Arm 2.
- **nf-core/rnaseq `-profile test`** run as an industry-standard fluency signal — ships tiny bundled data, runs on the laptop, no big compute needed.
- **Environments:** one **conda/mamba** env per arm to avoid dependency conflicts:
  - `env_bulk` — R/Bioconductor (DESeq2, edgeR, tximport, clusterProfiler, fgsea) + Salmon, STAR, subread, fastp, FastQC, MultiQC.
  - `env_singlecell` — Python/scverse (scanpy, scrublet, celltypist, harmonypy, scvi-tools) + Seurat via a small R env.
  - `env_longread` — minimap2, IsoQuant, Bambu (R), NanoPlot, seqtk; DTU R packages.
- **MultiQC** roll-ups where applicable.
- **README** with the pipeline DAG diagram, data provenance, and the explicit subsetting/downsampling disclosure.

---

## 10. Proposed repository layout

```
project2-rnaseq/
├── README.md
├── environments/
│   ├── env_bulk.yml
│   ├── env_singlecell.yml
│   └── env_longread.yml
├── config/
│   └── samples.tsv
├── data/                      # gitignored; populated by download scripts
├── scripts/
│   ├── download_sgnex.sh
│   ├── download_pbmc.sh
│   └── subsample.sh
├── arm1_bulk/
│   ├── Snakefile
│   ├── rules/
│   └── R/                     # DESeq2, edgeR, clusterProfiler, fgsea
├── arm2_singlecell/
│   ├── notebooks/             # scanpy_*.ipynb
│   └── seurat/                # seurat_parallel.R
├── arm3_longread/
│   ├── Snakefile
│   ├── rules/
│   └── R/                     # bambu, DTU
├── nfcore_demo/
│   └── run_test_profile.sh
└── results/                   # figures + tables per arm
```

---

## 11. Build order (milestones)

1. **Scaffold** — repo layout, env files, Snakemake skeletons, download/subsample scripts, README stub.
2. **Arm 1 bulk** — QC → Salmon → DESeq2/edgeR → functional; first figures.
3. **Arm 3 long-read** — minimap2 → IsoQuant/Bambu → DTU; complementarity figure with Arm 1.
4. **Arm 2 single-cell** — Scanpy notebook → annotation → Seurat parallel.
5. **nf-core test run** + README polish + DAG diagram.
6. **Final report** — figures assembled, short-vs-long story written up.

---

## 12. Known limitations & framing

- **Chr20–22 subset / downsampled reads:** disclosed as compute scoping; pipelines are genome- and depth-scalable.
- **STAR full-genome** not run locally (RAM); demonstrated on subset, noted as HPC-scalable.
- **Single-cell + long-read single-cell** not combined (niche, out of scope); noted as a possible extension.

---

## 13. Open items to confirm

- Exact SG-NEx cell-line pair (default: A549 vs K562).
- PBMC dataset choice vs a small tumor dataset for Arm 2.
- Whether direct-RNA ONT is included (adds RNA-modification angle) or cDNA only.
