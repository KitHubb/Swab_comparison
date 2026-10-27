# Current analysis workflow

Run these R Markdown files from the repository root in numeric order. Each file locates the repository root automatically and writes outputs below `result/final/`.

1. `00_Preprocessing.Rmd`
2. `01_Main_Supplement_Analysis.Rmd`
3. `02_Species_Composition.Rmd`
4. `03_Decontam_Sensitivity.Rmd` — full 55-sample contamination assessment,
   submitted Figure S1/S2 reproduction, and reviewer-requested residual-control profile.

The former Copan-forehead six-sample tutorial is archived at
`analysis/old/04_Decontam_Copan_Forehead_tutorial.Rmd`; it is not part of the
manuscript analysis.

`R/functions.R` contains shared plotting and taxonomic helper functions.
