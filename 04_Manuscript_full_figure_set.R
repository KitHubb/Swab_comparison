# Complete manuscript figure set
# Figure 1, Figure 2 (UpSet replacement), Figures S1-S3
# Existing analysis scripts are not modified.

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(ComplexUpset)
  library(pheatmap)
  library(RColorBrewer)
  library(ggpubr)
})

source("./Script/functions.R")

set.seed(42)

# -----------------------------------------------------------------------------
# 1. Data and common design
# -----------------------------------------------------------------------------

phy <- readRDS("./Phyloseq/phy_F270R240_260812_v4.rds")
phy <- subset_samples(phy, control_status == "sample") %>%
  prune_taxa(taxa_sums(.) > 0, .)

# Match the rarefied depth used in the existing analysis.
set.seed(42)
phy_rar <- rarefy_even_depth(
  phy,
  sample.size = min(sample_sums(phy)),
  rngseed = 42,
  replace = FALSE,
  verbose = FALSE
)

final_only <- identical(tolower(Sys.getenv("MANUSCRIPT_FINAL_ONLY", "false")), "true")
out_dir <- file.path(
  getwd(),
  if (final_only) "Final_main_supplement_figure_output" else "Manuscript_figure_set_output"
)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

swab_key <- data.frame(
  swab_type2 = c("Copan", "Puritan_gel", "Puritan_buffer", "Omni", "Eswab"),
  Type = factor(paste0("Type ", 1:5), levels = paste0("Type ", 1:5))
)

type_col <- setNames(brewer.pal(5, "Set2"), paste0("Type ", 1:5))
subject_col <- setNames(
  brewer.pal(5, "Set1"),
  sort(unique(as.character(sample_data(phy)$subject_id)))
)

theme_ms <- theme_bw(base_size = 9) +
  theme(
    axis.title = element_text(face = "bold", colour = "black"),
    axis.text = element_text(colour = "black"),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(5, 5, 5, 5)
  )

save_figure <- function(plot, stem, width, height) {
  # Apply consistent breathing room between manuscript panels while retaining
  # the requested overall figure dimensions.
  plot <- plot & theme(plot.margin = margin(8, 8, 8, 8, unit = "pt"))
  ggsave(file.path(out_dir, paste0(stem, ".png")), plot, width = width,
         height = height, dpi = 600, bg = "white", limitsize = FALSE)
  ggsave(file.path(out_dir, paste0(stem, ".tiff")), plot, width = width,
         height = height, dpi = 600, compression = "lzw", bg = "white",
         limitsize = FALSE)
}

meta_rar <- data.frame(sample_data(phy_rar), check.names = FALSE)
meta_rar$sample_name <- sample_names(phy_rar)
meta_rar <- left_join(meta_rar, swab_key, by = "swab_type2")

# -----------------------------------------------------------------------------
# 2. Figure 1: beta and alpha diversity across sampling systems
# -----------------------------------------------------------------------------

make_pcoa <- function(ps, site, method, label) {
  keep_samples <- sample_names(ps)[as.character(sample_data(ps)$site_group) == site]
  ps_site <- prune_samples(keep_samples, ps) %>%
    prune_taxa(taxa_sums(.) > 0, .)
  ord <- ordinate(ps_site, method = "PCoA", distance = method)
  eig <- ord$values$Relative_eig[1:2] * 100
  scores <- data.frame(ord$vectors[, 1:2], check.names = FALSE) %>%
    rownames_to_column("sample_name") %>%
    left_join(
      meta_rar %>% select(sample_name, subject_id, Type, swab_type2),
      by = "sample_name"
    )
  dist_obj <- phyloseq::distance(ps_site, method = method)
  md <- data.frame(sample_data(ps_site))
  perm <- adonis2(dist_obj ~ swab_type2, data = md, permutations = 9999)
  perm_r2 <- round(perm$R2[1], 3)
  perm_p <- perm$`Pr(>F)`[1]
  p_text <- ifelse(perm_p < 0.001, "<0.001", format(round(perm_p, 3), nsmall = 3))

  ggplot(scores, aes(Axis.1, Axis.2, colour = Type, shape = subject_id)) +
    geom_vline(xintercept = 0, colour = "grey80") +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_point(size = 2) +
    stat_ellipse(aes(group = Type), linewidth = 0.55, level = 0.95) +
    scale_colour_manual(values = type_col, drop = FALSE) +
    labs(
      x = NULL,
      y = paste0("PCoA 2 (", round(eig[2], 1), "%)"),
      subtitle = label,
      colour = "Sampling type",
      shape = "Participant"
    ) +
    annotate(
      "text", x = -Inf, y = -Inf, hjust = 0, vjust = -0.2, size = 2.7,
      label = paste0("\nPERMANOVA:\nR2 = ", perm_r2, ", p-value = ", p_text)
    ) +
    theme_ms +
    theme(plot.subtitle = element_text(size = 7, hjust = 0.5), aspect.ratio = 1)
}

alpha <- estimate_richness(phy_rar, measures = c("Shannon", "Observed")) %>%
  rownames_to_column("sample_name") %>%
  mutate(sample_name = sub("^X", "", sample_name)) %>%
  left_join(meta_rar, by = "sample_name")

make_sampling_alpha <- function(site) {
  alpha_site <- alpha %>%
    filter(site_group == site) %>%
    pivot_longer(c(Shannon, Observed), names_to = "Index", values_to = "Value") %>%
    mutate(Index = factor(Index, levels = c("Shannon", "Observed"),
                          labels = c("Shannon index", "Observed ASVs")))

  ggplot(alpha_site, aes(Type, Value, fill = Type)) +
    geom_boxplot(width = 0.65, outlier.shape = NA, colour = "black", linewidth = 0.45) +
    geom_point(
      aes(colour = Type, shape = subject_id),
      position = position_jitter(width = 0.08),
      size = 1.7
    ) +
    ggpubr::geom_pwc(
      method = "wilcox_test",
      label = "{p.adj.signif}",
      p.adjust.method = "BH",
      step.increase = 0.10,
      hide.ns = TRUE
    ) +
    facet_wrap(~Index, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = type_col, drop = FALSE) +
    scale_colour_manual(values = type_col, drop = FALSE) +
    labs(x = "Sampling type", y = NULL, title = site,
         fill = "Sampling type", colour = "Sampling type", shape = "Participant") +
    theme_ms +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(face = "bold"),
      strip.background = element_blank()
    )
}

f1_A <- make_pcoa(phy_rar, "Antecubital", "unifrac", "Unweighted UniFrac")
f1_B <- make_pcoa(phy_rar, "Antecubital", "wunifrac", "Weighted UniFrac")
f1_C <- make_pcoa(phy_rar, "Forehead", "unifrac", "Unweighted UniFrac")
f1_D <- make_pcoa(phy_rar, "Forehead", "wunifrac", "Weighted UniFrac")
f1_E <- make_sampling_alpha("Antecubital")
f1_F <- make_sampling_alpha("Forehead")

Figure1 <- (
  wrap_elements(full = (f1_A | f1_B)) /
  wrap_elements(full = (f1_C | f1_D)) /
  (f1_E | f1_F)
) +
  plot_layout(heights = c(1, 1, 0.62)) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D"))) &
  theme(plot.tag = element_text(size = 16, face = "bold"))
save_figure(Figure1, "Figure1_bacterial_diversity", 9, 10.2)

# Additional beta-diversity views following the existing beta_plot convention:
# participant = Set1 colour; swab type = point shape.
make_beta_subject_panel <- function(ps, site, method) {
  keep_samples <- sample_names(ps)[as.character(sample_data(ps)$site_group) == site]
  ps_site <- prune_samples(keep_samples, ps) %>% prune_taxa(taxa_sums(.) > 0, .)
  ord <- ordinate(ps_site, "PCoA", method)
  eig <- ord$values$Relative_eig[1:2] * 100
  scores <- data.frame(ord$vectors[, 1:2]) %>%
    rownames_to_column("sample_name") %>%
    left_join(meta_rar %>% select(sample_name, subject_id, swab_type2), by = "sample_name") %>%
    left_join(swab_key, by = "swab_type2")
  dist_obj <- phyloseq::distance(ps_site, method = method)
  md <- data.frame(sample_data(ps_site), check.names = FALSE)
  perm <- adonis2(dist_obj ~ subject_id, data = md, permutations = 9999)
  perm_r2 <- round(perm$R2[1], 3)
  perm_p <- perm$`Pr(>F)`[1]
  p_text <- ifelse(perm_p < 0.001, "<0.001", format(round(perm_p, 3), nsmall = 3))
  ggplot(scores, aes(Axis.1, Axis.2, colour = subject_id, shape = Type)) +
    geom_vline(xintercept = 0, colour = "grey85") +
    geom_hline(yintercept = 0, colour = "grey85") +
    geom_point(size = 2, alpha = 0.75) +
    stat_ellipse(aes(group = subject_id), linewidth = 0.55) +
    scale_colour_manual(values = subject_col) +
    labs(
      x = paste0("PCoA 1 (", round(eig[1], 1), "%)"),
      y = paste0("PCoA 2 (", round(eig[2], 1), "%)"),
      title = site,
      colour = "Participant",
      shape = "Swab type"
    ) +
    annotate(
      "text", x = -Inf, y = -Inf, hjust = 0, vjust = -0.2, size = 2.7,
      label = paste0("\nPERMANOVA:\nR2 = ", perm_r2, ", p-value = ", p_text)
    ) + theme_ms + theme(aspect.ratio = 1)
}

metric_labels <- c(
  bray = "Bray-Curtis", jaccard = "Jaccard",
  unifrac = "Unweighted UniFrac", wunifrac = "Weighted UniFrac"
)

if (!final_only) for (metric in names(metric_labels)) {
  beta_figure <-
    make_beta_subject_panel(phy_rar, "Antecubital", metric) |
    make_beta_subject_panel(phy_rar, "Forehead", metric) +
    plot_annotation(title = metric_labels[[metric]], tag_levels = "A") &
    theme(plot.tag = element_text(size = 14, face = "bold"))
  save_figure(beta_figure, paste0("Beta_PCoA_", metric), 7.5, 3.8)
}

# -----------------------------------------------------------------------------
# 3. Figure 2: bacterial community composition
# -----------------------------------------------------------------------------

phy_genus <- tax_glom(phy, taxrank = "Genus")
phy_genus_rel <- transform_sample_counts(phy_genus, function(x) x / sum(x))
major_genus_taxa <- taxa_names(phy_genus_rel)[
  taxa_sums(phy_genus_rel) / nsamples(phy_genus_rel) > 0.01
]

genus_tax_plot <- taxa_plot(
  melt = psmelt(phy_genus),
  taxa = "Genus",
  tax_otu = major_genus_taxa,
  x_axis = "swab_type2",
  phylum_or = "Actinobacteriota"
)

genus_composition_data <- genus_tax_plot$data$df %>%
  left_join(swab_key, by = "swab_type2")

make_genus_composition_panel <- function(site) {
  ggplot(
    filter(genus_composition_data, site_group == site),
    aes(Type, Abundance, fill = Genus)
  ) +
    geom_col(position = "fill", width = 0.88) +
    scale_fill_manual(values = genus_tax_plot$color) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
    labs(x = "Type", y = paste(site, "relative abundance (%)"), fill = "Genus") +
    guides(fill = guide_legend(keyheight = unit(4.8, "mm"), keywidth = unit(7, "mm"))) +
    theme_classic(base_size = 8) +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      legend.text = element_text(size = 6),
      legend.title = element_text(face = "bold"),
      legend.key.height = unit(4.8, "mm"),
      legend.key.width = unit(7, "mm")
    )
}

genus_rel_long <- psmelt(phy_genus_rel) %>%
  left_join(swab_key, by = "swab_type2") %>%
  mutate(Genus = as.character(Genus))

selected_skin_genera <- c("Cutibacterium", "Staphylococcus", "Streptococcus")

make_selected_genus_boxplot <- function(site) {
  ggplot(
    genus_rel_long %>%
      filter(site_group == site, Genus %in% selected_skin_genera),
    aes(Type, Abundance * 100, fill = Type)
  ) +
    geom_boxplot(width = 0.68, outlier.shape = NA, colour = "black") +
    geom_point(
      aes(shape = subject_id),
      colour = "black",
      position = position_jitter(width = 0.08),
      size = 1.5
    ) +
    facet_wrap(~Genus, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = type_col, drop = FALSE) +
    labs(x = "Type", y = paste(site, "relative abundance (%)"),
         fill = "Sampling type", shape = "Participant") +
    theme_classic(base_size = 8) +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "italic"),
      legend.position = "right"
    )
}

f2_A <- make_genus_composition_panel("Antecubital") +
  labs(tag = "A") + theme(plot.tag = element_text(size = 14, face = "bold"))
f2_B <- make_genus_composition_panel("Forehead") +
  labs(tag = "B") + theme(plot.tag = element_text(size = 14, face = "bold"))
f2_C <- make_selected_genus_boxplot("Antecubital") +
  labs(tag = "C") + theme(plot.tag = element_text(size = 14, face = "bold"))
f2_D <- make_selected_genus_boxplot("Forehead") +
  labs(tag = "D") + theme(plot.tag = element_text(size = 14, face = "bold"))

Figure2 <- (
  (f2_A | f2_B) /
  wrap_elements(full = f2_C) /
  wrap_elements(full = f2_D)
)

save_figure(Figure2, "Figure2_bacterial_community_composition", 8, 8.8)

if (!final_only) {
  write.csv(
    genus_composition_data,
    file.path(out_dir, "Figure2_genus_composition_data.csv"),
    row.names = FALSE
  )
  write.csv(
    genus_rel_long %>% filter(Genus %in% selected_skin_genera),
    file.path(out_dir, "Figure2_selected_genus_boxplot_data.csv"),
    row.names = FALSE
  )
}

# -----------------------------------------------------------------------------
# 4. Figure 3: UpSet plots at ASV, Genus, and Species levels
#    Taxa retained when mean relative abundance >=0.1% in each subset.
# -----------------------------------------------------------------------------

make_upset_data <- function(ps, site = NULL) {
  if (!is.null(site)) {
    keep_samples <- sample_names(ps)[as.character(sample_data(ps)$site_group) == site]
    ps <- prune_samples(keep_samples, ps) %>%
      prune_taxa(taxa_sums(.) > 0, .)
  }
  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
  keep <- taxa_names(ps_rel)[taxa_sums(ps_rel) / nsamples(ps_rel) >= 0.001]
  ps <- prune_taxa(keep, ps)
  otu <- data.frame(as(otu_table(ps), "matrix"), check.names = FALSE)
  if (!taxa_are_rows(ps)) otu <- t(otu)
  md <- data.frame(sample_data(ps), check.names = FALSE)
  md$sample_name <- sample_names(ps)
  present <- otu > 0
  result <- sapply(swab_key$swab_type2, function(sw) {
    ids <- rownames(md)[md$swab_type2 == sw]
    rowSums(present[, ids, drop = FALSE]) > 0
  })
  result <- data.frame(result, check.names = FALSE)
  colnames(result) <- as.character(swab_key$Type)
  result$taxon_id <- rownames(result)
  list(data = result, n_taxa = nrow(result))
}

make_upset <- function(obj, title) {
  ComplexUpset::upset(
    obj$data,
    intersect = paste0("Type ", 1:5),
    name = "Sampling system",
    min_size = 0,
    sort_sets = FALSE,
    sort_intersections_by = "cardinality",
    width_ratio = 0.2,
    base_annotations = list(
      "Intersection size" = intersection_size(counts = TRUE)
    )
  ) +
    labs(title = paste0(title, " (n=", obj$n_taxa, " taxa)")) +
    theme_bw(base_size = 8) +
    theme(plot.title = element_text(face = "bold"), panel.grid = element_blank())
}

make_upset_figure <- function(ps, level_name) {
  up_total <- make_upset_data(ps)
  up_ac <- make_upset_data(ps, "Antecubital")
  up_fh <- make_upset_data(ps, "Forehead")
  fig <- make_upset(up_total, "A  Total") /
    make_upset(up_ac, "B  Antecubital fossa") /
    make_upset(up_fh, "C  Forehead")
  save_figure(fig, paste0("Figure3_UpSet_", level_name), 7.5, 11)
  if (!final_only) {
    write.csv(up_total$data, file.path(out_dir, paste0("Figure3_", level_name, "_total.csv")), row.names = FALSE)
    write.csv(up_ac$data, file.path(out_dir, paste0("Figure3_", level_name, "_antecubital.csv")), row.names = FALSE)
    write.csv(up_fh$data, file.path(out_dir, paste0("Figure3_", level_name, "_forehead.csv")), row.names = FALSE)
  }
}

make_upset_figure(phy, "ASV")
make_upset_figure(tax_glom(phy, "Genus"), "Genus")
make_upset_figure(tax_glom(phy, "Species"), "Species")

# -----------------------------------------------------------------------------
# 5. Figure S1: participant-level composition at three taxonomic ranks
# -----------------------------------------------------------------------------

# Genus and Species use the existing taxa_plot() function directly so the
# phylum-specific palette and Other grouping exactly follow functions.R.
make_taxa_plot_composition <- function(ps, rank, mean_cutoff) {
  pg <- tax_glom(ps, taxrank = rank)
  pg_rel <- transform_sample_counts(pg, function(x) x / sum(x))
  keep <- taxa_names(pg_rel)[taxa_sums(pg_rel) / nsamples(pg_rel) > mean_cutoff]
  tx <- taxa_plot(
    melt = psmelt(pg),
    taxa = rank,
    tax_otu = keep,
    x_axis = "swab_type2",
    phylum_or = "Actinobacteriota"
  )
  plot_data <- tx$data$df %>% left_join(swab_key, by = "swab_type2")
  ggplot(plot_data, aes(Type, Abundance, fill = .data[[rank]])) +
    geom_col(position = "fill", width = 0.92) +
    facet_grid(site_group ~ subject_id) +
    scale_fill_manual(values = tx$color, drop = FALSE) +
    guides(fill = guide_legend(keyheight = unit(4.4, "mm"), keywidth = unit(7, "mm"))) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
    labs(x = NULL, y = "Relative abundance (%)", fill = rank) +
    theme_classic(base_size = 7) +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(size = 6),
      strip.background = element_rect(fill = "white", colour = "black"),
      strip.text = element_text(face = "bold"),
      legend.text = element_text(size = 6),
      legend.key.height = unit(4.4, "mm"),
      legend.key.width = unit(7, "mm")
    )
}

# taxa_plot() cannot take Phylum as both its grouping and taxon column.
# These colours reproduce its F3.generate_colors() mapping for one taxon per
# phylum (Reds/Blues/Purples/Greens and grey for Other).
make_phylum_composition <- function(ps, mean_cutoff = 0.01) {
  pg <- tax_glom(ps, "Phylum")
  pg_rel <- transform_sample_counts(pg, function(x) x / sum(x))
  keep <- taxa_names(pg_rel)[taxa_sums(pg_rel) / nsamples(pg_rel) > mean_cutoff]
  long <- psmelt(pg)
  keep_phyla <- unique(as.character(tax_table(pg)[keep, "Phylum"]))
  long$Phylum <- ifelse(as.character(long$Phylum) %in% keep_phyla,
                        as.character(long$Phylum), "Other")
  long <- long %>% left_join(swab_key, by = "swab_type2")
  phylum_palette_name <- c(
    Actinobacteriota = "Reds", Actinomycetota = "Reds",
    Firmicutes = "Blues", Bacillota = "Blues",
    Bacteroidota = "Purples", Proteobacteria = "Greens",
    Pseudomonadota = "Greens", Fusobacteriota = "YlOrBr"
  )
  phyla <- setdiff(sort(unique(long$Phylum)), "Other")
  cols <- setNames(vapply(phyla, function(x) {
    pal <- if (x %in% names(phylum_palette_name)) phylum_palette_name[[x]] else "BrBG"
    rev(brewer.pal(9, pal))[5]
  }, character(1)), phyla)
  cols <- c(cols, Other = "#D3D3D3")
  long$Phylum <- factor(long$Phylum, levels = names(cols))
  ggplot(long, aes(Type, Abundance, fill = Phylum)) +
    geom_col(position = "fill", width = 0.92) +
    facet_grid(site_group ~ subject_id) +
    scale_fill_manual(values = cols, drop = FALSE) +
    guides(fill = guide_legend(keyheight = unit(4.4, "mm"), keywidth = unit(7, "mm"))) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
    labs(x = NULL, y = "Relative abundance (%)", fill = "Phylum") +
    theme_classic(base_size = 7) +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(size = 6),
      strip.background = element_rect(fill = "white", colour = "black"),
      strip.text = element_text(face = "bold"),
      legend.text = element_text(size = 6),
      legend.key.height = unit(4.4, "mm"),
      legend.key.width = unit(7, "mm")
    )
}

s1_A <- make_phylum_composition(phy, 0.01)
s1_B <- make_taxa_plot_composition(phy, "Genus", 0.01)
s1_C <- make_taxa_plot_composition(phy, "Species", 0.01)
FigureS1 <- (s1_A / s1_B / s1_C) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))
save_figure(FigureS1, "FigureS1_participant_composition", 10.7, 10.8)

# -----------------------------------------------------------------------------
# 6. Figure S2: distance heatmaps by site
# -----------------------------------------------------------------------------

subject_annotation_col <- subject_col
swab_annotation_col <- type_col

make_heatmap <- function(ps, site, method) {
  keep_samples <- sample_names(ps)[as.character(sample_data(ps)$site_group) == site]
  ps_site <- prune_samples(keep_samples, ps) %>% prune_taxa(taxa_sums(.) > 0, .)
  mat <- as.matrix(phyloseq::distance(ps_site, method = method))
  md <- data.frame(sample_data(ps_site), check.names = FALSE) %>%
    rownames_to_column("sample_name") %>% left_join(swab_key, by = "swab_type2")
  ann <- data.frame(Subject = md$subject_id, Swab_type = md$Type)
  rownames(ann) <- md$sample_name
  pheatmap(
    mat, silent = TRUE, clustering_method = "complete",
    annotation_col = ann, annotation_row = ann,
    annotation_colors = list(Subject = subject_annotation_col, Swab_type = swab_annotation_col),
    show_rownames = TRUE, show_colnames = TRUE, fontsize = 6,
    border_color = NA,
    main = paste(site, metric_labels[[method]])
  )$gtable
}

s2_A <- wrap_elements(full = make_heatmap(phy_rar, "Antecubital", "unifrac"))
s2_B <- wrap_elements(full = make_heatmap(phy_rar, "Antecubital", "wunifrac"))
s2_C <- wrap_elements(full = make_heatmap(phy_rar, "Forehead", "unifrac"))
s2_D <- wrap_elements(full = make_heatmap(phy_rar, "Forehead", "wunifrac"))
FigureS2 <- (s2_A | s2_B) / (s2_C | s2_D) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))
save_figure(FigureS2, "FigureS2_UniFrac_heatmaps", 10.5, 8.5)

# Save Bray-Curtis, Jaccard, unweighted UniFrac and weighted UniFrac separately.
if (!final_only) for (metric in names(metric_labels)) {
  heatmap_figure <-
    wrap_elements(full = make_heatmap(phy_rar, "Antecubital", metric)) |
    wrap_elements(full = make_heatmap(phy_rar, "Forehead", metric)) +
    plot_annotation(title = metric_labels[[metric]], tag_levels = "A") &
    theme(plot.tag = element_text(size = 14, face = "bold"))
  save_figure(
    heatmap_figure,
    paste0("FigureS2_distance_clustering_", metric),
    7.2,
    3.7
  )
}

# -----------------------------------------------------------------------------
# 7. Figure S3: participant-level alpha and beta diversity
# -----------------------------------------------------------------------------

alpha_s3 <- estimate_richness(phy_rar, measures = c("Shannon", "Observed")) %>%
  rownames_to_column("sample_name") %>%
  mutate(sample_name = sub("^X", "", sample_name)) %>%
  left_join(meta_rar, by = "sample_name") %>%
  pivot_longer(c(Shannon, Observed), names_to = "Index", values_to = "Value") %>%
  mutate(Index = factor(Index, levels = c("Shannon", "Observed"),
                        labels = c("Shannon index", "Observed ASVs")))

make_participant_alpha <- function(site) {
  ggplot(
    filter(alpha_s3, site_group == site),
    aes(subject_id, Value, fill = subject_id)
  ) +
    geom_boxplot(outlier.shape = NA, width = 0.62, colour = "black", linewidth = 0.45) +
    geom_point(
      aes(shape = Type),
      colour = "black",
      position = position_jitter(width = 0.10),
      size = 1.5
    ) +
    ggpubr::geom_pwc(
      method = "wilcox_test",
      label = "{p.adj.signif}",
      p.adjust.method = "BH",
      step.increase = 0.10,
      hide.ns = TRUE
    ) +
    facet_wrap(~Index, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = subject_col, drop = FALSE) +
    labs(x = "Participant", y = paste(site, "diversity"),
         fill = "Participant", shape = "Sampling type") +
    theme_classic(base_size = 8) +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "right",
      legend.key.height = unit(3, "mm")
    )
}

s3_A <- make_participant_alpha("Antecubital")
s3_B <- make_participant_alpha("Forehead")

alpha_s3_statistics <- ggpubr::compare_means(
  Value ~ subject_id,
  data = alpha_s3,
  group.by = c("site_group", "Index"),
  method = "wilcox.test",
  p.adjust.method = "BH"
)

if (!final_only) {
  write.csv(
    alpha_s3_statistics,
    file.path(out_dir, "FigureS3_alpha_pairwise_Wilcoxon_BH.csv"),
    row.names = FALSE
  )
}

make_subject_pcoa <- function(site, method) {
  keep_samples <- sample_names(phy_rar)[
    as.character(sample_data(phy_rar)$site_group) == site
  ]
  ps <- prune_samples(keep_samples, phy_rar) %>% prune_taxa(taxa_sums(.) > 0, .)
  ord <- ordinate(ps, "PCoA", method)
  eig <- ord$values$Relative_eig[1:2] * 100
  sc <- data.frame(ord$vectors[, 1:2]) %>% rownames_to_column("sample_name") %>%
    left_join(meta_rar %>% select(sample_name, subject_id, Type), by = "sample_name")
  dist_obj <- phyloseq::distance(ps, method = method)
  md <- data.frame(sample_data(ps), check.names = FALSE)
  perm <- adonis2(dist_obj ~ subject_id, data = md, permutations = 9999)
  perm_r2 <- round(perm$R2[1], 3)
  perm_p <- perm$`Pr(>F)`[1]
  p_text <- ifelse(perm_p < 0.001, "<0.001", format(round(perm_p, 3), nsmall = 3))
  p <- ggplot(sc, aes(Axis.1, Axis.2, colour = subject_id, shape = Type)) +
    geom_vline(xintercept = 0, colour = "grey80") +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_point(size = 1.8) +
    stat_ellipse(aes(group = subject_id), linewidth = 0.5) +
    scale_colour_manual(values = subject_col) +
    labs(x = paste0("PCoA 1 (", round(eig[1], 1), "%)"),
         y = paste0("PCoA 2 (", round(eig[2], 1), "%)"),
         subtitle = ifelse(method == "wunifrac", "Weighted UniFrac", "Unweighted UniFrac"),
         colour = "Participant", shape = "Swab type") +
    annotate(
      "text", x = -Inf, y = -Inf, hjust = 0, vjust = -0.2, size = 2.5,
      label = paste0("\nPERMANOVA:\nR2 = ", perm_r2, ", p-value = ", p_text)
    ) + theme_ms + theme(aspect.ratio = 1)

  if (method == "unifrac") {
    p <- p + theme(legend.position = "none")
  }
  p
}

s3_C <- make_subject_pcoa("Antecubital", "unifrac") | make_subject_pcoa("Antecubital", "wunifrac")
s3_D <- make_subject_pcoa("Forehead", "unifrac") | make_subject_pcoa("Forehead", "wunifrac")
FigureS3 <- (s3_A | s3_B) /
  wrap_elements(full = s3_C) /
  wrap_elements(full = s3_D) +
  plot_layout(heights = c(0.72, 1, 1)) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D"))) &
  theme(plot.tag = element_text(size = 16, face = "bold"))
save_figure(FigureS3, "FigureS3_participant_alpha_beta", 11.5, 9.5)

if (!final_only) {
  write.csv(swab_key, file.path(out_dir, "Swab_type_key.csv"), row.names = FALSE)
}
cat("\nManuscript figure output folder:\n", normalizePath(out_dir), "\n")
print(list.files(out_dir, full.names = TRUE))
