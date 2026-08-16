# Final manuscript figures only
# Outputs Figure 1-3 and Figure S1-S3 as PNG/TIFF without auxiliary plots/CSVs.

Sys.setenv(MANUSCRIPT_FINAL_ONLY = "true")
source(file.path(getwd(), "04_Manuscript_full_figure_set.R"), chdir = TRUE)
