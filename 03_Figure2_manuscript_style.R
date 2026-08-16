# Figure generation matched to the visual format of the CoST manuscript
# Existing analysis scripts are not modified by this file.

library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)


# -----------------------------------------------------------------------------
# 1. Input and output
# -----------------------------------------------------------------------------

phy <- readRDS("./Phyloseq/phy_F270R240_260812_v4.rds")
phy <- subset_samples(phy, control_status == "sample") %>%
  prune_taxa(taxa_sums(.) > 0, .)

figure_output_dir <- file.path(getwd(), "Figure2_manuscript_style_output")

if (!dir.exists(figure_output_dir)) {
  dir.create(figure_output_dir, recursive = TRUE)
}


# Mapping follows the order used throughout the existing analysis.
swab_key <- data.frame(
  swab_type2 = c(
    "Copan", "Eswab", "Omni", "Puritan_buffer", "Puritan_gel"
  ),
  Type = factor(paste0("Type ", 1:5), levels = paste0("Type ", 1:5))
)

write.csv(
  swab_key,
  file.path(figure_output_dir, "Swab_type_key.csv"),
  row.names = FALSE
)


type_colours <- c(
  "Type 1" = "#F8766D",
  "Type 2" = "#A3A500",
  "Type 3" = "#00BF7D",
  "Type 4" = "#00B0F6",
  "Type 5" = "#E76BF3"
)


theme_manuscript <- theme_classic(base_size = 10, base_family = "sans") +
  theme(
    axis.title = element_text(face = "bold", colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.45),
    axis.ticks = element_line(colour = "black", linewidth = 0.4),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(7, 7, 7, 7)
  )


# -----------------------------------------------------------------------------
# 2. Swab performance dataset
# -----------------------------------------------------------------------------

otu_sample <- data.frame(as(otu_table(phy), "matrix"), check.names = FALSE)

if (taxa_are_rows(phy)) {
  otu_sample <- t(otu_sample) %>%
    data.frame(check.names = FALSE)
}

meta_sample <- data.frame(sample_data(phy), check.names = FALSE)
meta_sample$sample_name <- sample_names(phy)

performance_data <- data.frame(
  sample_name = rownames(otu_sample),
  total_reads = rowSums(otu_sample),
  observed_ASV = rowSums(otu_sample > 0),
  check.names = FALSE
) %>%
  left_join(
    meta_sample %>%
      select(sample_name, subject_id, site_group, swab_type2),
    by = "sample_name"
  ) %>%
  left_join(swab_key, by = "swab_type2")


otu_long <- otu_sample %>%
  rownames_to_column("sample_name") %>%
  pivot_longer(
    cols = -sample_name,
    names_to = "taxon_id",
    values_to = "abundance"
  ) %>%
  left_join(
    meta_sample %>%
      select(sample_name, subject_id, site_group, swab_type2),
    by = "sample_name"
  ) %>%
  mutate(detected = abundance > 0) %>%
  group_by(subject_id, site_group) %>%
  mutate(union_ASV = n_distinct(taxon_id[detected])) %>%
  group_by(subject_id, site_group, swab_type2, sample_name, union_ASV) %>%
  summarise(
    detected_ASV = sum(detected),
    recovery_rate = detected_ASV / first(union_ASV),
    .groups = "drop"
  )


performance_data <- performance_data %>%
  left_join(
    otu_long %>%
      select(sample_name, union_ASV, recovery_rate),
    by = "sample_name"
  )


# -----------------------------------------------------------------------------
# 3. Statistical summaries
# -----------------------------------------------------------------------------

performance_summary <- performance_data %>%
  group_by(site_group, Type, swab_type2) %>%
  summarise(
    n = n(),
    median_reads = median(total_reads),
    IQR_reads = IQR(total_reads),
    median_observed_ASV = median(observed_ASV),
    IQR_observed_ASV = IQR(observed_ASV),
    median_recovery_rate = median(recovery_rate),
    IQR_recovery_rate = IQR(recovery_rate),
    .groups = "drop"
  )


run_friedman <- function(data, outcome) {
  bind_rows(lapply(split(data, data$site_group), function(x) {
    fit <- friedman.test(x[[outcome]], x$Type, x$subject_id)
    data.frame(
      site_group = x$site_group[1],
      outcome = outcome,
      statistic = unname(fit$statistic),
      df = unname(fit$parameter),
      p_value = fit$p.value
    )
  }))
}


friedman_results <- bind_rows(
  run_friedman(performance_data, "total_reads"),
  run_friedman(performance_data, "observed_ASV"),
  run_friedman(performance_data, "recovery_rate")
) %>%
  mutate(p_adjust_BH = p.adjust(p_value, method = "BH"))


write.csv(
  performance_data,
  file.path(figure_output_dir, "Swab_performance_raw.csv"),
  row.names = FALSE
)

write.csv(
  performance_summary,
  file.path(figure_output_dir, "Swab_performance_summary.csv"),
  row.names = FALSE
)

write.csv(
  friedman_results,
  file.path(figure_output_dir, "Swab_performance_Friedman_test.csv"),
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# 4. Main Figure: manuscript-style 2 x 3 panel layout
# -----------------------------------------------------------------------------

make_performance_panel <- function(data, site, outcome, y_label) {
  plot_data <- data %>% filter(site_group == site)
  
  ggplot(
    plot_data,
    aes(x = Type, y = .data[[outcome]], group = subject_id)
  ) +
    geom_line(colour = "grey65", linewidth = 0.45, alpha = 0.8) +
    geom_point(
      aes(fill = Type),
      shape = 21,
      size = 2.8,
      stroke = 0.55,
      colour = "black"
    ) +
    stat_summary(
      aes(group = Type),
      fun = median,
      geom = "crossbar",
      width = 0.52,
      linewidth = 0.55,
      colour = "black"
    ) +
    scale_fill_manual(values = type_colours, drop = FALSE) +
    labs(x = "Sampling system", y = y_label) +
    theme_manuscript +
    theme(legend.position = "none")
}


plot_A <- make_performance_panel(
  performance_data, "Antecubital", "total_reads", "Antecubital\nTotal reads"
)
plot_B <- make_performance_panel(
  performance_data, "Antecubital", "observed_ASV", "Antecubital\nObserved ASVs"
)
plot_C <- make_performance_panel(
  performance_data, "Antecubital", "recovery_rate", "Antecubital\nTaxa recovery rate"
) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

plot_D <- make_performance_panel(
  performance_data, "Forehead", "total_reads", "Forehead\nTotal reads"
)
plot_E <- make_performance_panel(
  performance_data, "Forehead", "observed_ASV", "Forehead\nObserved ASVs"
)
plot_F <- make_performance_panel(
  performance_data, "Forehead", "recovery_rate", "Forehead\nTaxa recovery rate"
) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))


figure_swab_performance <-
  (plot_A | plot_B | plot_C) /
  (plot_D | plot_E | plot_F) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(
      size = 18,
      face = "bold",
      colour = "black"
    )
  )


ggsave(
  file.path(figure_output_dir, "Figure_swab_performance_manuscript_style.png"),
  figure_swab_performance,
  width = 12,
  height = 7.8,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(figure_output_dir, "Figure_swab_performance_manuscript_style.tiff"),
  figure_swab_performance,
  width = 12,
  height = 7.8,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# -----------------------------------------------------------------------------
# 5. Major-genus detection prevalence: manuscript-style A-B panels
# -----------------------------------------------------------------------------

phy_genus <- tax_glom(phy, taxrank = "Genus") %>%
  prune_taxa(taxa_sums(.) > 0, .)

phy_genus_rel <- transform_sample_counts(phy_genus, function(x) x / sum(x))
major_genus_taxa <- names(sort(taxa_sums(phy_genus_rel), decreasing = TRUE))[1:10]

otu_genus <- data.frame(
  as(otu_table(phy_genus), "matrix"),
  check.names = FALSE
)

if (taxa_are_rows(phy_genus)) {
  otu_genus <- t(otu_genus) %>%
    data.frame(check.names = FALSE)
}

tax_genus <- data.frame(tax_table(phy_genus), check.names = FALSE) %>%
  rownames_to_column("taxon_id") %>%
  mutate(
    Genus = ifelse(
      is.na(Genus) | Genus == "",
      "Unclassified",
      as.character(Genus)
    )
  )

genus_prevalence <- otu_genus %>%
  rownames_to_column("sample_name") %>%
  pivot_longer(
    cols = -sample_name,
    names_to = "taxon_id",
    values_to = "abundance"
  ) %>%
  filter(taxon_id %in% major_genus_taxa) %>%
  left_join(tax_genus %>% select(taxon_id, Genus), by = "taxon_id") %>%
  left_join(
    meta_sample %>%
      select(sample_name, site_group, swab_type2),
    by = "sample_name"
  ) %>%
  left_join(swab_key, by = "swab_type2") %>%
  group_by(site_group, Type, swab_type2, Genus) %>%
  summarise(
    detection_prevalence = mean(abundance > 0),
    .groups = "drop"
  )


genus_order <- tax_genus$Genus[
  match(major_genus_taxa, tax_genus$taxon_id)
]

make_detection_panel <- function(data, site) {
  data %>%
    filter(site_group == site) %>%
    mutate(Genus = factor(Genus, levels = rev(genus_order))) %>%
    ggplot(aes(x = Type, y = Genus, fill = detection_prevalence)) +
    geom_tile(colour = "white", linewidth = 0.7) +
    geom_text(
      aes(label = scales::percent(detection_prevalence, accuracy = 1)),
      size = 3,
      colour = "black"
    ) +
    scale_fill_gradient(
      low = "white",
      high = "#2C7FB8",
      limits = c(0, 1),
      labels = scales::percent_format(accuracy = 1)
    ) +
    labs(x = "Sampling system", y = NULL, fill = "Detection\nprevalence") +
    theme_manuscript +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.45)
    )
}


figure_genus_detection <- (
  make_detection_panel(genus_prevalence, "Antecubital") +
    labs(title = "Antecubital fossa") |
  make_detection_panel(genus_prevalence, "Forehead") +
    labs(title = "Forehead")
) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(
      size = 18,
      face = "bold",
      colour = "black"
    ),
    plot.title = element_text(face = "bold", hjust = 0)
  )


ggsave(
  file.path(figure_output_dir, "Figure_major_genus_detection_manuscript_style.png"),
  figure_genus_detection,
  width = 10,
  height = 5.8,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(figure_output_dir, "Figure_major_genus_detection_manuscript_style.tiff"),
  figure_genus_detection,
  width = 10,
  height = 5.8,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

write.csv(
  genus_prevalence,
  file.path(figure_output_dir, "Major_genus_detection_prevalence.csv"),
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# 6. Console output
# -----------------------------------------------------------------------------

cat("\nFigure output folder:\n")
cat(normalizePath(figure_output_dir), "\n\n")
cat("Friedman test results:\n")
print(friedman_results)
cat("\nGenerated files:\n")
print(list.files(figure_output_dir, full.names = TRUE))
