# Analysis-ready data

`phyloseq_analysis.rds` is the current post-decontamination object used by the main analysis. When preprocessing is rerun, the initial imported object is written as `phyloseq_imported.rds`. `phyloseq_predecontam.rds` and `concentration_metadata.csv` support contamination sensitivity checks. `taxonomy_blast.qza` supports the putative species-level supplementary analysis. The TSV files provide portable abundance, taxonomy, metadata, and quality-control summaries.

Earlier processed objects are retained in `old/`, where their original locations are recorded.
