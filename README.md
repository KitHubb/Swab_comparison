# Comparative Evaluation of Skin Microbiome Sampling Systems

## Introduction

This repository contains the R project, analysis code, processed data objects, tables, and figures used to reproduce the analyses reported in the manuscript *“Comparative Evaluation of Skin Microbiome Sampling Systems.”* The study compares five integrated swab–collection-medium systems for V1–V3 16S rRNA gene profiling of the antecubital fossa and forehead in a repeated-measures design.

## Preprocessing

Raw paired-end reads were processed with the Nextflow workflow [amplicon_16S_v1v3_qiime_nf](https://github.com/KitHubb/amplicon_16S_v1v3_qiime_nf).

```bash
conda activate nextflow_nf

cd <ANALYSIS_PROJECT_DIRECTORY>

nextflow run \
  <PATH_TO_PIPELINE>/main.nf \
  -profile singularity \
  -params-file <PATH_TO_Q20_PARAMETER_FILE> \
  --reads '<PATH_TO_READS>/*_{1,2}.fastq.gz' \
  --run_label <RUN_LABEL> \
  --outdir <PATH_TO_OUTPUT_DIRECTORY> \
  --classifier <PATH_TO_SILVA_CLASSIFIER> \
  --taxonomy_label SILVA \
  --metadata <PATH_TO_QIIME2_METADATA> \
  --trimm_optimal true \
  --trimm_combinations <PATH_TO_TRUNCATION_PARAMETER_FILE> \
  --diversity_enabled false \
  -work-dir <PATH_TO_WORK_DIRECTORY> \
  -resume
```

### Decontamination

Prevalence-based contaminant identification was evaluated across decontam thresholds by comparing the numbers and proportions of ASVs and reads retained or removed from biological samples and negative controls. The supporting sensitivity workflow is distributed as the R package [decontamSensitivity](https://github.com/KitHubb/decontamSensitivity).

## Downstream analysis

QIIME 2 outputs were converted into `phyloseq` objects and analysed in R. The workflow includes singleton removal, rarefaction, alpha and beta diversity, taxonomic composition, overlap among sampling systems, and PERMANOVA. Detailed analysis and figure-generation code is provided in the files below.

## Core files

| File | Purpose |
|---|---|
| `01_Preprocessing.Rmd` | Metadata integration, decontamination, quality checks, and construction of analysis-ready `phyloseq` objects |
| `02_Main_Supplement_Figure_Tables_script.Rmd` | Statistical analyses and generation of the manuscript's main and supplementary figures and tables |
| `03_supple_figuresS3.rmd` | Species level results |
| `Figures/` | Final main and supplementary figures |
| `Tables/` | Final tables and numerical results |
| `Phyloseq/` | Processed `phyloseq` objects used by the R workflow |
| `metadata/` | Analysis metadata |

## Data availability

The sequence data from this study have been submitted to the NCBI BioProject under accession number [PRJNA1075916](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1075916).

## Citation

Citation information will be added after publication of the associated article.
