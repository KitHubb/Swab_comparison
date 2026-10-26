# Real-data decontam execution check

Archived bundle. Original repository path: `examples/decontam_realdata_check/`.

This folder checks whether the prevalence, frequency, and combined methods from
`decontam::isContaminant()` run successfully on the supplied DT_Swab phyloseq
object. It does not use functions from `decontamSensitivity` and is not wired
into a package build or an automated Git workflow.

## Inputs

- `../../Phyloseq/phy_F270R240_260812_v3.rds`
- `../../../Paper/HN00182797_Original_Sample_QC_and_NGS_Library_QC.xlsx`
- Excel sheet: `Sequencing_Metadata`
- Frequency concentration: `Original_QC_Conc_ng_uL`
- Negative-control indicator: phyloseq sample-data column `is.neg`, or the five
  missing `control_status` rows when `is.neg` is absent in v3

The Excel sample names are required to match the phyloseq sample names exactly.
The original-QC concentration is used because it is positive for all 55 samples.
The post-PCR library concentration is not used because it contains two zero
values and measures a different stage of the workflow.

## Run

From this folder:

```powershell
python prepare_concentration_csv.py
Rscript run_decontam_realdata_check.R
Rscript render_threshold_comparison_report.R
```

The scripts validate sample matching, control classification, positive
concentrations, detailed-score/direct-classification agreement, and thresholds
from 0.1 through 0.7 in 0.1 steps. Outputs are written under `output/`.
The report command writes a self-contained
`output/decontam_threshold_comparison_v3.html` with threshold, read-retention,
agreement, flag-pattern, and taxonomy views.

## Interpretation limitation

The analysis uses the pre-decontam `v3.rds` object so prevalence, frequency,
and combined methods are compared from the same unfiltered feature table.

The original-QC concentration range is also narrow (0.40-0.56 ng/uL; CV about
8.2%). A successful frequency run confirms software compatibility and input
alignment, but does not by itself establish strong frequency-based separation.
