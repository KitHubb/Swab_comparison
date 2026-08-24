# Publication-oriented validation of the canonical singleton analysis.
# Starts from the post-decontam object and does not overwrite prior outputs.

project_library <- file.path(getwd(), "R_library")
if (dir.exists(project_library)) .libPaths(c(project_library, .libPaths()))

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

set.seed(42)
output_dir <- file.path(getwd(), "Results", "canonical_singleton_20260824")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

phy_post_decontam <- readRDS("Phyloseq/phy_F270R240_260812_v4.rds")
phy_true <- subset_samples(phy_post_decontam, control_status == "sample")
phy_true <- prune_taxa(taxa_sums(phy_true) > 0, phy_true)
phy_singleton <- prune_taxa(taxa_sums(phy_true) > 1, phy_true)

filter_summary <- tibble(
  samples = nsamples(phy_singleton),
  asvs_before = ntaxa(phy_true),
  asvs_after = ntaxa(phy_singleton),
  asvs_removed = ntaxa(phy_true) - ntaxa(phy_singleton),
  reads_before = sum(sample_sums(phy_true)),
  reads_after = sum(sample_sums(phy_singleton)),
  reads_removed = reads_before - reads_after,
  reads_retained_pct = 100 * reads_after / reads_before,
  rarefaction_depth = min(sample_sums(phy_singleton))
)
write.csv(filter_summary, file.path(output_dir, "Table_filter_and_sample_summary.csv"), row.names = FALSE)

set.seed(42)
phy_rar <- rarefy_even_depth(
  phy_singleton,
  sample.size = min(sample_sums(phy_singleton)),
  rngseed = 42,
  replace = FALSE,
  verbose = FALSE
)
metadata <- data.frame(sample_data(phy_rar), check.names = FALSE)

# Hypotheses: each variable explains no community variation for each distance.
# Participant is tested with unrestricted permutations; repeated-measure factors
# are tested using permutations restricted within participant.
distance_methods <- c(
  weighted_unifrac = "wunifrac",
  unweighted_unifrac = "unifrac",
  bray_curtis = "bray"
)
model_variables <- c("subject_id", "site_group", "swab_type2")

run_permanova <- function(distance_name, distance_method, variable) {
  distance_object <- phyloseq::distance(phy_rar, method = distance_method)
  model_formula <- as.formula(paste("distance_object ~", variable))
  blocked <- variable != "subject_id"
  set.seed(42)
  fit <- if (blocked) {
    adonis2(model_formula, data = metadata, permutations = 9999,
            strata = metadata$subject_id)
  } else {
    adonis2(model_formula, data = metadata, permutations = 9999)
  }
  tibble(
    distance = distance_name,
    explanatory_variable = variable,
    response = "community dissimilarity",
    n_samples = nrow(metadata),
    n_participants = n_distinct(metadata$subject_id),
    df = fit$Df[1],
    pseudo_F = fit$F[1],
    R2 = fit$R2[1],
    p_value = fit$`Pr(>F)`[1],
    permutations = 9999L,
    permutation_scheme = ifelse(blocked, "restricted within participant", "unrestricted")
  )
}

permanova_results <- bind_rows(lapply(names(distance_methods), function(distance_name) {
  bind_rows(lapply(model_variables, function(variable) {
    run_permanova(distance_name, distance_methods[[distance_name]], variable)
  }))
}))
write.csv(permanova_results, file.path(output_dir, "Table_PERMANOVA_univariable_all_distances.csv"), row.names = FALSE)

# PERMANOVA assumes comparable multivariate dispersion. This diagnostic tests
# dispersion among sampling systems with permutations restricted by participant.
run_dispersion <- function(distance_name, distance_method) {
  distance_object <- phyloseq::distance(phy_rar, method = distance_method)
  dispersion <- betadisper(distance_object, group = metadata$swab_type2,
                           type = "median", bias.adjust = TRUE)
  permutation_design <- permute::how(nperm = 9999, blocks = metadata$subject_id)
  set.seed(42)
  fit <- permutest(dispersion, permutations = permutation_design)
  tibble(
    distance = distance_name,
    explanatory_variable = "swab_type2",
    n_samples = nrow(metadata),
    n_groups = n_distinct(metadata$swab_type2),
    df = fit$tab$Df[1],
    F = fit$tab$F[1],
    p_value = fit$tab$`Pr(>F)`[1],
    permutations = 9999L,
    permutation_scheme = "restricted within participant"
  )
}

dispersion_results <- bind_rows(lapply(names(distance_methods), function(distance_name) {
  run_dispersion(distance_name, distance_methods[[distance_name]])
}))
write.csv(dispersion_results, file.path(output_dir, "Table_PERMDISP_sampling_system.csv"), row.names = FALSE)

# Paired alpha-diversity contrasts. Effect direction is type_2 minus type_1.
alpha <- estimate_richness(phy_rar, measures = c("Shannon", "Observed")) |>
  rownames_to_column("sample_name") |>
  mutate(sample_name = sub("^X", "", sample_name)) |>
  left_join(
    metadata |>
      rownames_to_column("sample_name") |>
      mutate(sample_name = sub("^X", "", sample_name)),
    by = "sample_name"
  )
sampling_types <- sort(unique(as.character(alpha$swab_type2)))
sampling_pairs <- combn(sampling_types, 2, simplify = FALSE)

paired_alpha_test <- function(site, index, pair) {
  paired <- alpha |>
    filter(site_group == site) |>
    select(subject_id, swab_type2, all_of(index)) |>
    pivot_wider(names_from = swab_type2, values_from = all_of(index)) |>
    filter(!is.na(.data[[pair[1]]]), !is.na(.data[[pair[2]]]))
  difference <- paired[[pair[2]]] - paired[[pair[1]]]
  ranks <- rank(abs(difference[difference != 0]))
  signs <- sign(difference[difference != 0])
  rank_biserial <- if (length(ranks)) sum(ranks * signs) / sum(ranks) else 0
  fit <- suppressWarnings(wilcox.test(
    paired[[pair[2]]], paired[[pair[1]]], paired = TRUE, exact = FALSE,
    conf.int = TRUE, conf.level = 0.95
  ))
  tibble(
    site = site,
    response = index,
    type_1 = pair[1],
    type_2 = pair[2],
    n_pairs = nrow(paired),
    median_type_1 = median(paired[[pair[1]]]),
    median_type_2 = median(paired[[pair[2]]]),
    median_paired_difference = median(difference),
    hodges_lehmann_difference = unname(fit$estimate),
    ci_low = fit$conf.int[1],
    ci_high = fit$conf.int[2],
    matched_rank_biserial = rank_biserial,
    p_value = fit$p.value
  )
}

alpha_results <- bind_rows(lapply(unique(as.character(alpha$site_group)), function(site) {
  bind_rows(lapply(c("Shannon", "Observed"), function(index) {
    bind_rows(lapply(sampling_pairs, function(pair) paired_alpha_test(site, index, pair)))
  }))
})) |>
  group_by(site, response) |>
  mutate(
    p_adjusted_bonferroni = p.adjust(p_value, method = "bonferroni"),
    p_adjusted_BH = p.adjust(p_value, method = "BH")
  ) |>
  ungroup()
write.csv(alpha_results, file.path(output_dir, "Table_alpha_pairwise_paired_effect_sizes.csv"), row.names = FALSE)

sanity_checks <- tibble(
  check = c(
    "all true samples retained",
    "minimum rarefaction depth is 5876",
    "all rarefied samples have equal depth",
    "all alpha contrasts have five participant pairs"
  ),
  passed = c(
    nsamples(phy_singleton) == 50,
    min(sample_sums(phy_singleton)) == 5876,
    length(unique(sample_sums(phy_rar))) == 1,
    all(alpha_results$n_pairs == 5)
  )
)
write.csv(sanity_checks, file.path(output_dir, "validation_sanity_checks.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

if (!all(sanity_checks$passed)) stop("One or more validation sanity checks failed.")
