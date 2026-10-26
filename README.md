# Comparative Evaluation of Skin Microbiome Sampling Systems

## Introduction

This repository contains the R project, analysis code, analysis-ready data objects, tables, and figures used to reproduce the manuscript analyses. The study compares five integrated swab–collection-medium systems for V1–V3 16S rRNA gene profiling of the antecubital fossa and forehead in a repeated-measures design.

## Repository structure

| Path | Purpose |
|---|---|
| `server/` | Nextflow configuration and server-side execution notes |
| `data/process/` | Analysis-ready metadata, phyloseq objects, and exported intermediate tables |
| `data/raw/` | Raw or external pipeline inputs; intentionally excluded from Git |
| `analysis/final/` | Current reproducible preprocessing and manuscript analysis Rmd files |
| `analysis/old/` | Earlier and exploratory analyses retained with their original paths |
| `result/final/` | Current manuscript figures, tables, and rendered reports |
| `result/old/` | Earlier results retained on GitHub with their original paths documented |
| `docs/` | Project documentation and manuscript-location notes |

Every `old/` directory contains an `ORIGINAL_PATHS.tsv` and a README so an earlier analysis can be related to the paths used when it was originally run.

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

### Blast modification 
Putative species-level supplementation using blast (using the "species_missing_rescue option"") [qiime2_blast_nf](https://github.com/KitHubb/qiime2_blast_nf).

```
nextflow run /data/software/nextflow/QIIME_blast/main.nf \
  --backend native \
  --repseq_qza <PATH_TO_INPUT_DIRECTORY>/05_dada2/rep-seqs.qza \
  --taxonomy_qza  <PATH_TO_INPUT_DIRECTORY>/06_taxonomy/taxonomy_SILVA.qza \
  --outdir  <PATH_TO_INPUT_DIRECTORY>/09_BLAST/taxonomy_SILVA_blast.qza \
  --blast_db /data/Reference/BLAST/16S/20241203/16S_ribosomal_RNA \
  --taxdump_dir /data/Reference/BLAST/taxdump \
  --taxonomy_profile generic \
  --blast_taxonomy_sif <BLAST_sif> \
  -work-dir <PATH_TO_WORK_DIRECTORY> \
  -resume

```


### Decontamination

Prevalence-based contaminant identification was evaluated across decontam thresholds by comparing the numbers and proportions of ASVs and reads retained or removed from biological samples and negative controls. The supporting sensitivity workflow is distributed as the R package [decontamSensitivity](https://github.com/KitHubb/decontamSensitivity).

## Analysis order

Run the current analyses from the repository root in this order:

1. `analysis/final/00_Preprocessing.Rmd`
2. `analysis/final/01_Main_Supplement_Analysis.Rmd`
3. `analysis/final/02_Species_Composition.Rmd`
4. `analysis/final/03_Decontam_Sensitivity.Rmd`

The downstream workflow starts from `data/process/phyloseq_analysis.rds`. It includes singleton removal, rarefaction, repeated-measures statistics, alpha and beta diversity, taxonomic composition, overlap among sampling systems, and PERMANOVA.

## Current core files

| File | Purpose |
|---|---|
| `analysis/final/00_Preprocessing.Rmd` | Metadata integration, decontamination, quality checks, and construction of analysis-ready objects |
| `analysis/final/01_Main_Supplement_Analysis.Rmd` | Repeated-measures analyses and generation of main and supplementary figures and tables |
| `analysis/final/02_Species_Composition.Rmd` | Putative species-level supplementary composition analysis |
| `analysis/final/03_Decontam_Sensitivity.Rmd` | Decontam threshold and control-profile sensitivity report |
| `result/final/figures/` | Current main and supplementary figures |
| `result/final/tables/` | Current numerical results and audit tables |
| `data/process/` | Processed inputs used by the R workflow |

## Data availability

The sequence data from this study have been submitted to the NCBI BioProject under accession number [PRJNA1075916](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1075916).

## Citation

Citation information will be added after publication of the associated article.
