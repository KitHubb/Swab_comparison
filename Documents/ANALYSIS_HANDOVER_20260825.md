# DT Swab sensitivity-analysis handover

- Date: 2026-08-25
- Base commit: `dc9851c`
- Status: canonical singleton validation completed; decontam-threshold benchmark blocked by pre-decontam control metadata

## Completed validation

The canonical post-decontam threshold 0.5 plus singleton condition retained 50 samples,
1,780 ASVs, and 847,074 reads. It removed 75 singleton ASVs and 75 reads from the
post-decontam true-sample object and retained 99.991% of reads. The fixed rarefaction
depth was 5,876 reads.

Sampling-system PERMDISP was not significant for weighted UniFrac (p=0.9151),
unweighted UniFrac (p=0.9418), or Bray-Curtis (p=0.7993). Thus, there was no evidence
that the non-significant sampling-system PERMANOVA results were driven by unequal
multivariate dispersion under the canonical singleton condition.

Univariable PERMANOVA effect sizes under the singleton condition were:

| Distance | Participant R2 | Site R2 | Sampling-system R2 | Sampling-system p |
|---|---:|---:|---:|---:|
| Weighted UniFrac | 0.4343 | 0.2742 | 0.0244 | 0.8027 |
| Unweighted UniFrac | 0.1985 | 0.0530 | 0.0615 | 0.8826 |
| Bray-Curtis | 0.4080 | 0.1155 | 0.0300 | 0.9947 |

The validation code is `07_Canonical_singleton_statistical_validation.R`.

## Planned decontam-threshold benchmark

`08_Decontam_filter_benchmark.Rmd` implements the prespecified 27-condition design:

- prevalence decontam thresholds 0.1 through 0.9;
- no additional rare-feature filtering;
- singleton removal;
- removal of ASVs with total count <=10.

For each condition it is designed to export data retention, participant/site/system
PERMANOVA, sampling-system PERMDISP, composition similarity, shared ASV/Genus read
fraction, and system-exclusive read fraction. Each decontam threshold is reconstructed
by rerunning `isContaminant()` as in `01_Preprocessing.Rmd`; the code does not manually
threshold the returned p-values.

## Blocking data issue

The benchmark must not be interpreted or reported yet. Its reconstruction sanity
check intentionally failed at threshold 0.5:

- reconstructed threshold 0.5/no-filter: 1,956 ASVs and 924,394 true-sample reads;
- saved v4 threshold 0.5 object: 1,855 ASVs and 847,149 true-sample reads.

Inspection showed that the five negative controls in
`Phyloseq/phy_F270R240_260812_v3.rds` have `control_status = NA`, whereas true samples
are labelled `sample`. Consequently, directly setting
`is.neg = control_status == "control"` produces 50 FALSE and 5 NA values, and decontam
does not reconstruct the original contaminant calls. The saved v4 object contains the
correct `control_status`, `is.neg`, `swab_type2`, and `SampleID` metadata.

Before rerunning the benchmark, restore the v3 control metadata from the original
metadata-generation step or `metadata/DT_metadata.tsv`. Do not infer control status
solely from an ad hoc sample-name rule without checking the preprocessing record.
After restoration, threshold 0.5/no-filter must reproduce exactly 1,855 ASVs,
847,149 reads, and 50 true samples before the remaining 26 combinations are accepted.

## Existing filtering-sensitivity interpretation

The existing threshold 0.5 x nine-filter analysis remains suitable for supplementary
reporting. No-filter, singleton, count <=2, and count <=10 conditions preserve the
sampling-system conclusion while retaining all 50 samples. Strong prevalence and
mean-abundance filters alter sample inclusion at the fixed rarefaction depth and must
therefore be presented as exploratory boundary conditions rather than evidence of
improved concordance.

## Git and environment

`R_library/` is a local package library and is excluded in `.gitignore`. It must not be
committed. Generated incomplete benchmark output should also remain uncommitted until
the reconstruction sanity checks pass and the report renders successfully.
