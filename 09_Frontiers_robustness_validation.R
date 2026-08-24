# Frontiers-oriented robustness validation for the canonical analysis.
# This script is independent of the existing manuscript Rmd files.
# All permutation-based tests use 9,999 permutations and seed 42.

project_library <- file.path(getwd(), "R_library")
if (dir.exists(project_library)) .libPaths(c(project_library, .libPaths()))

suppressPackageStartupMessages({
  library(phyloseq)
  library(ape)
  library(vegan)
  library(permute)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
})

set.seed(42)
permutations_n <- 9999L
target_depth <- 5876L
out_dir <- file.path(getwd(), "Frontiers_robustness_output")
fig_dir <- file.path(out_dir, "Figures")
tab_dir <- file.path(out_dir, "Tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

otu_matrix <- function(ps) {
  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu <- t(otu)
  otu
}

prepare_unifrac_tree <- function(ps) {
  tree <- phy_tree(ps)
  polytomy_present <- !ape::is.binary(tree)
  edges_before <- nrow(tree$edge)
  if (polytomy_present) {
    tree <- ape::multi2di(tree, random = FALSE)
    phy_tree(ps) <- tree
  }
  stopifnot(ape::is.binary(phy_tree(ps)))
  attr(ps, "tree_audit") <- list(
    polytomy_present = polytomy_present,
    edges_before = edges_before,
    edges_after = nrow(phy_tree(ps)$edge),
    added_zero_length_edges = nrow(phy_tree(ps)$edge) - edges_before
  )
  ps
}

phy_post <- readRDS("Phyloseq/phy_F270R240_260812_v4.rds")
phy_true <- subset_samples(phy_post, control_status == "sample")
phy_true <- prune_taxa(taxa_sums(phy_true) > 0, phy_true)
phy_singleton <- prune_taxa(taxa_sums(phy_true) > 1, phy_true)
stopifnot(
  nsamples(phy_singleton) == 50L,
  ntaxa(phy_singleton) == 1780L,
  sum(sample_sums(phy_singleton)) == 847074,
  min(sample_sums(phy_singleton)) == target_depth
)

set.seed(42)
phy_rar <- rarefy_even_depth(
  phy_singleton,
  sample.size = target_depth,
  rngseed = 42,
  replace = FALSE,
  trimOTUs = TRUE,
  verbose = FALSE
)
metadata <- data.frame(sample_data(phy_rar), check.names = FALSE)

distance_methods <- c(
  weighted_unifrac = "wunifrac",
  unweighted_unifrac = "unifrac",
  bray_curtis = "bray"
)
distance_objects <- lapply(distance_methods, function(method) {
  work <- if (method %in% c("wunifrac", "unifrac")) prepare_unifrac_tree(phy_rar) else phy_rar
  phyloseq::distance(work, method = method)
})

# -------------------------------------------------------------------------
# 1. Leave-one-participant-out PERMANOVA and sampling-system PERMDISP
# -------------------------------------------------------------------------

run_loo_permanova <- function(ps, distance_name, method, variable, omitted) {
  if (method %in% c("wunifrac", "unifrac")) ps <- prepare_unifrac_tree(ps)
  dist_obj <- phyloseq::distance(ps, method = method)
  md <- data.frame(sample_data(ps), check.names = FALSE)
  blocked <- variable != "subject_id"
  model_formula <- as.formula(paste("dist_obj ~", variable))
  set.seed(42)
  fit <- if (blocked) {
    adonis2(model_formula, data = md, permutations = permutations_n,
            strata = md$subject_id)
  } else {
    adonis2(model_formula, data = md, permutations = permutations_n)
  }
  tibble(
    omitted_participant = omitted,
    distance = distance_name,
    variable = variable,
    n_samples = nrow(md),
    n_participants = n_distinct(md$subject_id),
    df = fit$Df[1],
    pseudo_F = fit$F[1],
    R2 = fit$R2[1],
    p_value = fit$`Pr(>F)`[1],
    permutations = permutations_n,
    permutation_scheme = ifelse(blocked, "restricted within participant", "unrestricted")
  )
}

run_loo_permdisp <- function(ps, distance_name, method, omitted) {
  if (method %in% c("wunifrac", "unifrac")) ps <- prepare_unifrac_tree(ps)
  dist_obj <- phyloseq::distance(ps, method = method)
  md <- data.frame(sample_data(ps), check.names = FALSE)
  dispersion <- betadisper(dist_obj, md$swab_type2, type = "median", bias.adjust = TRUE)
  design <- permute::how(nperm = permutations_n, blocks = md$subject_id)
  set.seed(42)
  fit <- permutest(dispersion, permutations = design)
  tibble(
    omitted_participant = omitted,
    distance = distance_name,
    variable = "swab_type2",
    n_samples = nrow(md),
    F = fit$tab$F[1],
    p_value = fit$tab$`Pr(>F)`[1],
    permutations = permutations_n,
    permutation_scheme = "restricted within participant"
  )
}

participants <- sort(unique(as.character(metadata$subject_id)))
loo_permanova <- bind_rows(lapply(participants, function(omitted) {
  keep_samples <- rownames(metadata)[as.character(metadata$subject_id) != omitted]
  ps <- prune_samples(keep_samples, phy_rar)
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  bind_rows(lapply(names(distance_methods), function(distance_name) {
    bind_rows(lapply(c("subject_id", "site_group", "swab_type2"), function(variable) {
      run_loo_permanova(ps, distance_name, distance_methods[[distance_name]], variable, omitted)
    }))
  }))
}))

loo_permdisp <- bind_rows(lapply(participants, function(omitted) {
  keep_samples <- rownames(metadata)[as.character(metadata$subject_id) != omitted]
  ps <- prune_samples(keep_samples, phy_rar)
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  bind_rows(lapply(names(distance_methods), function(distance_name) {
    run_loo_permdisp(ps, distance_name, distance_methods[[distance_name]], omitted)
  }))
}))

loo_summary <- loo_permanova %>%
  group_by(distance, variable) %>%
  summarise(
    analyses = n(),
    R2_min = min(R2),
    R2_max = max(R2),
    p_min = min(p_value),
    p_max = max(p_value),
    significant_n = sum(p_value < 0.05),
    .groups = "drop"
  )

write.csv(loo_permanova, file.path(tab_dir, "Table1_leave_one_participant_out_PERMANOVA.csv"), row.names = FALSE)
write.csv(loo_permdisp, file.path(tab_dir, "Table2_leave_one_participant_out_PERMDISP.csv"), row.names = FALSE)
write.csv(loo_summary, file.path(tab_dir, "Table3_leave_one_participant_out_summary.csv"), row.names = FALSE)

# -------------------------------------------------------------------------
# 2. Matched participant x site distance benchmark
# -------------------------------------------------------------------------

pair_table <- function(dist_obj, distance_name) {
  mat <- as.matrix(dist_obj)
  idx <- which(upper.tri(mat), arr.ind = TRUE)
  out <- tibble(
    sample_1 = rownames(mat)[idx[, 1]],
    sample_2 = colnames(mat)[idx[, 2]],
    dissimilarity = mat[idx],
    distance = distance_name
  )
  md <- metadata %>% rownames_to_column("sample_name")
  out %>%
    left_join(md %>% select(sample_name, subject_id, site_group, swab_type2),
              by = c("sample_1" = "sample_name")) %>%
    rename(subject_1 = subject_id, site_1 = site_group, system_1 = swab_type2) %>%
    left_join(md %>% select(sample_name, subject_id, site_group, swab_type2),
              by = c("sample_2" = "sample_name")) %>%
    rename(subject_2 = subject_id, site_2 = site_group, system_2 = swab_type2) %>%
    mutate(
      comparison = case_when(
        subject_1 == subject_2 & site_1 == site_2 & system_1 != system_2 ~
          "Within participant-site: different systems",
        subject_1 == subject_2 & site_1 != site_2 & system_1 == system_2 ~
          "Within participant: different sites, same system",
        subject_1 != subject_2 & site_1 == site_2 & system_1 == system_2 ~
          "Between participants: same site and system",
        TRUE ~ "Other"
      )
    ) %>%
    filter(comparison != "Other")
}

matched_distances <- bind_rows(lapply(names(distance_objects), function(distance_name) {
  pair_table(distance_objects[[distance_name]], distance_name)
}))

matched_distance_summary <- matched_distances %>%
  group_by(distance, comparison) %>%
  summarise(
    pair_count = n(),
    median = median(dissimilarity),
    q1 = quantile(dissimilarity, 0.25),
    q3 = quantile(dissimilarity, 0.75),
    mean = mean(dissimilarity),
    sd = sd(dissimilarity),
    .groups = "drop"
  )

within_block_summary <- matched_distances %>%
  filter(comparison == "Within participant-site: different systems") %>%
  mutate(block = paste(subject_1, site_1, sep = "__")) %>%
  group_by(distance, block, subject_id = subject_1, site_group = site_1) %>%
  summarise(median_within_block_distance = median(dissimilarity), .groups = "drop")

write.csv(matched_distances, file.path(tab_dir, "Table4_matched_pairwise_distances.csv"), row.names = FALSE)
write.csv(matched_distance_summary, file.path(tab_dir, "Table5_matched_distance_summary.csv"), row.names = FALSE)
write.csv(within_block_summary, file.path(tab_dir, "Table6_within_block_distance_by_participant_site.csv"), row.names = FALSE)

# Pairwise distances are not independent; this figure and table are descriptive.
comparison_labels <- c(
  "Within participant-site: different systems" = "Within block\nSystems differ",
  "Within participant: different sites, same system" = "Between sites\nSame participant",
  "Between participants: same site and system" = "Between participants\nSame site"
)
distance_labels <- c(
  bray_curtis = "Bray-Curtis",
  unweighted_unifrac = "Unweighted UniFrac",
  weighted_unifrac = "Weighted UniFrac"
)
p_distance <- matched_distances %>%
  mutate(plot_comparison = factor(
    unname(comparison_labels[comparison]),
    levels = unname(comparison_labels)
  )) %>%
  ggplot(aes(plot_comparison, dissimilarity, fill = plot_comparison)) +
  geom_violin(scale = "width", alpha = 0.6, colour = NA) +
  geom_boxplot(width = 0.16, outlier.shape = NA, colour = "black") +
  facet_wrap(~distance, scales = "free_y", labeller = as_labeller(distance_labels)) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = NULL, y = "Dissimilarity") +
  theme_bw(base_size = 9) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 8, lineheight = 0.95)
  )

ggsave(file.path(fig_dir, "Figure1_matched_distance_benchmark.png"),
       p_distance, width = 10, height = 5.5, dpi = 600)
ggsave(file.path(fig_dir, "Figure1_matched_distance_benchmark.tiff"),
       p_distance, width = 10, height = 5.5, dpi = 600, compression = "lzw")

# -------------------------------------------------------------------------
# 3. Participant x site block-specific ASV and Genus agreement
# -------------------------------------------------------------------------

prepare_rank_object <- function(ps, rank_name) {
  if (rank_name == "ASV") return(ps)
  tax <- as(tax_table(ps), "matrix")
  genus <- as.character(tax[, "Genus"])
  genus[is.na(genus) | genus == ""] <- "Unclassified_genus"
  tax[, "Genus"] <- genus
  tax_table(ps) <- tax_table(tax)
  tax_glom(ps, taxrank = "Genus", NArm = FALSE)
}

block_agreement <- function(ps, rank_name) {
  work <- prepare_rank_object(ps, rank_name)
  otu <- otu_matrix(work)
  md <- data.frame(sample_data(work), check.names = FALSE)
  blocks <- md %>%
    rownames_to_column("sample_name") %>%
    distinct(subject_id, site_group)

  bind_rows(lapply(seq_len(nrow(blocks)), function(i) {
    block_md <- md[md$subject_id == blocks$subject_id[i] &
                     md$site_group == blocks$site_group[i], , drop = FALSE]
    samples <- rownames(block_md)
    block_otu <- otu[, samples, drop = FALSE]
    present <- block_otu > 0
    shared <- rownames(block_otu)[rowSums(present) == ncol(block_otu)]
    union_feature_ids <- rownames(block_otu)[rowSums(present) > 0]
    intersection_over_union_value <- length(shared) / length(union_feature_ids)
    tibble(
      taxonomic_level = rank_name,
      subject_id = blocks$subject_id[i],
      site_group = blocks$site_group[i],
      systems_present = n_distinct(block_md$swab_type2),
      union_features = length(union_feature_ids),
      shared_all_five_features = length(shared),
      intersection_over_union = intersection_over_union_value,
      shared_all_five_read_pct = 100 * sum(block_otu[shared, , drop = FALSE]) / sum(block_otu)
    )
  }))
}

block_agreement_results <- bind_rows(
  block_agreement(phy_singleton, "ASV"),
  block_agreement(phy_singleton, "Genus")
)

block_agreement_summary <- block_agreement_results %>%
  group_by(taxonomic_level) %>%
  summarise(
    blocks = n(),
    median_shared_features = median(shared_all_five_features),
    shared_features_min = min(shared_all_five_features),
    shared_features_max = max(shared_all_five_features),
    median_intersection_over_union = median(intersection_over_union),
    median_shared_read_pct = median(shared_all_five_read_pct),
    shared_read_pct_min = min(shared_all_five_read_pct),
    shared_read_pct_max = max(shared_all_five_read_pct),
    .groups = "drop"
  )

write.csv(block_agreement_results, file.path(tab_dir, "Table7_participant_site_block_agreement.csv"), row.names = FALSE)
write.csv(block_agreement_summary, file.path(tab_dir, "Table8_participant_site_block_agreement_summary.csv"), row.names = FALSE)

p_shared <- block_agreement_results %>%
  ggplot(aes(site_group, shared_all_five_read_pct, colour = subject_id,
             group = subject_id)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 2) +
  facet_wrap(~taxonomic_level) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Anatomical site", y = "Reads assigned to features shared by all five systems (%)",
       colour = "Participant") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "Figure2_block_specific_shared_read_fraction.png"),
       p_shared, width = 8, height = 4.5, dpi = 600)
ggsave(file.path(fig_dir, "Figure2_block_specific_shared_read_fraction.tiff"),
       p_shared, width = 8, height = 4.5, dpi = 600, compression = "lzw")

# -------------------------------------------------------------------------
# 4. Repeated-measures omnibus tests for alpha diversity and read depth
# -------------------------------------------------------------------------

alpha <- estimate_richness(phy_rar, measures = c("Shannon", "Observed")) %>%
  rownames_to_column("sample_name") %>%
  mutate(sample_name = sub("^X", "", sample_name)) %>%
  left_join(
    metadata %>%
      rownames_to_column("sample_name") %>%
      mutate(sample_name = sub("^X", "", sample_name)),
    by = "sample_name"
  )

run_friedman <- function(data, response_name, site) {
  subset_data <- data %>% filter(site_group == site)
  fit <- friedman.test(
    y = subset_data[[response_name]],
    groups = subset_data$swab_type2,
    blocks = subset_data$subject_id
  )
  tibble(
    site_group = site,
    response = response_name,
    n_participants = n_distinct(subset_data$subject_id),
    n_systems = n_distinct(subset_data$swab_type2),
    statistic = unname(fit$statistic),
    df = unname(fit$parameter),
    p_value = fit$p.value,
    design = "Friedman repeated-measures omnibus test"
  )
}

alpha_friedman <- bind_rows(lapply(unique(as.character(alpha$site_group)), function(site) {
  bind_rows(
    run_friedman(alpha, "Shannon", site),
    run_friedman(alpha, "Observed", site)
  )
}))

read_depth <- tibble(
  sample_name = sample_names(phy_singleton),
  post_filter_reads = as.numeric(sample_sums(phy_singleton))
) %>%
  mutate(sample_name = sub("^X", "", sample_name)) %>%
  left_join(
    data.frame(sample_data(phy_singleton), check.names = FALSE) %>%
      rownames_to_column("sample_name") %>%
      mutate(sample_name = sub("^X", "", sample_name)),
    by = "sample_name"
  )

read_depth_friedman <- bind_rows(lapply(unique(as.character(read_depth$site_group)), function(site) {
  run_friedman(read_depth, "post_filter_reads", site)
}))

write.csv(alpha_friedman, file.path(tab_dir, "Table9_alpha_diversity_Friedman_tests.csv"), row.names = FALSE)
write.csv(read_depth_friedman, file.path(tab_dir, "Table10_postfilter_reads_Friedman_tests.csv"), row.names = FALSE)

# -------------------------------------------------------------------------
# 5. Reproducibility and validation
# -------------------------------------------------------------------------

tree_before <- phy_tree(phy_rar)
tree_resolved_ps <- prepare_unifrac_tree(phy_rar)
tree_after <- phy_tree(tree_resolved_ps)
tree_audit <- tibble(
  tips = ape::Ntip(tree_before),
  internal_nodes_before = tree_before$Nnode,
  internal_nodes_after = tree_after$Nnode,
  edges_before = nrow(tree_before$edge),
  edges_after = nrow(tree_after$edge),
  zero_length_edges_before = sum(tree_before$edge.length == 0),
  zero_length_edges_after = sum(tree_after$edge.length == 0),
  rooted_before = ape::is.rooted(tree_before),
  binary_before = ape::is.binary(tree_before),
  binary_after = ape::is.binary(tree_after),
  tip_set_preserved = setequal(tree_before$tip.label, tree_after$tip.label)
)
write.csv(tree_audit, file.path(tab_dir, "Table11_UniFrac_tree_audit.csv"), row.names = FALSE)

sanity_checks <- tibble(
  check = c(
    "canonical singleton sample count",
    "canonical singleton ASV count",
    "canonical singleton read count",
    "rarefaction depth",
    "five leave-one-participant-out runs",
    "all block analyses contain five systems",
    "ten participant-site blocks per taxonomic level",
    "all permutation tests use 9999",
    "UniFrac tree is binary after zero-length resolution"
  ),
  passed = c(
    nsamples(phy_singleton) == 50L,
    ntaxa(phy_singleton) == 1780L,
    sum(sample_sums(phy_singleton)) == 847074,
    length(unique(sample_sums(phy_rar))) == 1L && unique(sample_sums(phy_rar)) == target_depth,
    n_distinct(loo_permanova$omitted_participant) == 5L,
    all(block_agreement_results$systems_present == 5L),
    all(table(block_agreement_results$taxonomic_level) == 10L),
    all(loo_permanova$permutations == permutations_n) &&
      all(loo_permdisp$permutations == permutations_n),
    tree_audit$binary_after && tree_audit$tip_set_preserved
  )
)
write.csv(sanity_checks, file.path(tab_dir, "validation_sanity_checks.csv"), row.names = FALSE)
if (!all(sanity_checks$passed)) stop("One or more Frontiers robustness sanity checks failed")

writeLines(c(
  "Frontiers robustness analysis",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "Input: Phyloseq/phy_F270R240_260812_v4.rds",
  "Canonical condition: decontam prevalence threshold 0.5 plus singleton removal",
  "Rarefaction depth: 5876",
  "Permutation tests: 9999 permutations",
  "Random seed: 42",
  "Matched-pair distance comparisons are descriptive because pairwise distances are non-independent."
), file.path(out_dir, "README_results.txt"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))

message("Frontiers robustness validation completed successfully.")
