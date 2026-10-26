script_directory <- function() {
  command <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command, value = TRUE)
  if (length(file_arg)) {
    return(dirname(sub("^--file=", "", file_arg[1L])))
  }
  getwd()
}

analysis_dir <- script_directory()
source(file.path(analysis_dir, "run_decontam_realdata_check.R"), local = TRUE)

required_packages <- c("ggplot2", "scales", "base64enc")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing report packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

ggplot2::theme_set(
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(color = "#52616f"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      strip.text = ggplot2::element_text(face = "bold")
    )
)

method_levels <- c("prevalence", "frequency", "combined")
method_labels <- c(
  prevalence = "Prevalence",
  frequency = "Frequency",
  combined = "Combined"
)
method_colors <- c(
  prevalence = "#007c91",
  frequency = "#d97706",
  combined = "#6d5bd0"
)
threshold_summary$method <- factor(
  threshold_summary$method,
  levels = method_levels
)

save_plot_uri <- function(plot, width, height) {
  plot_path <- tempfile(fileext = ".png")
  on.exit(unlink(plot_path), add = TRUE)
  ggplot2::ggsave(
    plot_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 150,
    bg = "white"
  )
  base64enc::dataURI(file = plot_path, mime = "image/png")
}

html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

html_table <- function(data, digits = 3L, class_name = "data-table") {
  format_cell <- function(x) {
    if (is.numeric(x)) {
      out <- ifelse(
        is.na(x),
        "",
        format(round(x, digits), trim = TRUE, scientific = FALSE, big.mark = ",")
      )
    } else if (is.logical(x)) {
      out <- ifelse(is.na(x), "", ifelse(x, "Yes", "No"))
    } else {
      out <- ifelse(is.na(x), "", as.character(x))
    }
    html_escape(out)
  }
  header <- paste0(
    "<thead><tr>",
    paste0("<th>", html_escape(names(data)), "</th>", collapse = ""),
    "</tr></thead>"
  )
  body_rows <- vapply(seq_len(nrow(data)), function(i) {
    paste0(
      "<tr>",
      paste0(
        "<td>",
        vapply(data[i, , drop = FALSE], format_cell, character(1)),
        "</td>",
        collapse = ""
      ),
      "</tr>"
    )
  }, character(1))
  paste0(
    '<div class="table-wrap"><table class="', class_name, '">',
    header,
    "<tbody>",
    paste(body_rows, collapse = ""),
    "</tbody></table></div>"
  )
}

plot_flagged <- ggplot2::ggplot(
  threshold_summary,
  ggplot2::aes(
    x = threshold,
    y = contaminant_features,
    color = method,
    shape = method
  )
) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(size = 2.4) +
  ggplot2::scale_color_manual(
    values = method_colors,
    labels = method_labels
  ) +
  ggplot2::scale_shape_manual(
    values = c(prevalence = 16, frequency = 17, combined = 15),
    labels = method_labels
  ) +
  ggplot2::scale_x_continuous(breaks = thresholds) +
  ggplot2::scale_y_continuous(labels = scales::comma) +
  ggplot2::labs(
    title = "Flagged ASVs across thresholds",
    subtitle = "v3 pre-decontamination feature table; 2,108 ASVs",
    x = "decontam threshold",
    y = "Flagged ASVs",
    color = NULL,
    shape = NULL
  )

retention_long <- rbind(
  data.frame(
    method = threshold_summary$method,
    threshold = threshold_summary$threshold,
    sample_group = "All reads",
    retained_pct = 100 - threshold_summary$reads_removed_total_pct
  ),
  data.frame(
    method = threshold_summary$method,
    threshold = threshold_summary$threshold,
    sample_group = "Biological samples",
    retained_pct = 100 - threshold_summary$reads_removed_biological_pct
  ),
  data.frame(
    method = threshold_summary$method,
    threshold = threshold_summary$threshold,
    sample_group = "Negative controls",
    retained_pct = 100 - threshold_summary$reads_removed_control_pct
  )
)
retention_long$sample_group <- factor(
  retention_long$sample_group,
  levels = c("All reads", "Biological samples", "Negative controls")
)
plot_retained <- ggplot2::ggplot(
  retention_long,
  ggplot2::aes(x = threshold, y = retained_pct, color = method, shape = method)
) +
  ggplot2::geom_hline(yintercept = 90, color = "#bcc6cf", linewidth = 0.4) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.0) +
  ggplot2::facet_wrap(~sample_group, nrow = 1) +
  ggplot2::scale_color_manual(values = method_colors, labels = method_labels) +
  ggplot2::scale_shape_manual(
    values = c(prevalence = 16, frequency = 17, combined = 15),
    labels = method_labels
  ) +
  ggplot2::scale_x_continuous(breaks = thresholds) +
  ggplot2::scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%")
  ) +
  ggplot2::labs(
    title = "Reads retained after removing flagged ASVs",
    subtitle = "Sharp drops reveal thresholds that remove a large fraction of the library",
    x = "decontam threshold",
    y = "Reads retained",
    color = NULL,
    shape = NULL
  )

overlap_summary$comparison <- paste0(
  method_labels[overlap_summary$method_a],
  " vs ",
  method_labels[overlap_summary$method_b]
)
plot_overlap <- ggplot2::ggplot(
  overlap_summary,
  ggplot2::aes(x = threshold, y = jaccard, color = comparison, shape = comparison)
) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::scale_x_continuous(breaks = thresholds) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = scales::percent_format(accuracy = 1)
  ) +
  ggplot2::scale_color_manual(values = c("#007c91", "#6d5bd0", "#d97706")) +
  ggplot2::labs(
    title = "Agreement between flagged-ASV sets",
    subtitle = "Jaccard similarity: shared flags divided by the union of flags",
    x = "decontam threshold",
    y = "Jaccard similarity",
    color = NULL,
    shape = NULL
  )

taxonomy_df <- if (is.null(taxonomy)) {
  data.frame(feature_id = feature_ids, Genus = "Unavailable")
} else {
  data.frame(
    feature_id = rownames(taxonomy),
    Genus = as.character(taxonomy[["Genus"]]),
    stringsAsFactors = FALSE
  )
}
taxonomy_df$Genus[is.na(taxonomy_df$Genus) | taxonomy_df$Genus == ""] <-
  "Unclassified"
feature_read_totals <- rowSums(counts)

flagged_taxonomy <- do.call(rbind, lapply(method_levels, function(method) {
  do.call(rbind, lapply(seq_along(thresholds), function(i) {
    flags <- method_results[[method]]$flags[, i]
    selected_ids <- feature_ids[flags]
    if (!length(selected_ids)) return(NULL)
    selected_genera <- taxonomy_df$Genus[
      match(selected_ids, taxonomy_df$feature_id)
    ]
    aggregate(
      cbind(flagged_reads = feature_read_totals[selected_ids], asv_count = 1) ~
        Genus,
      data = data.frame(
        Genus = selected_genera,
        flagged_reads = unname(feature_read_totals[selected_ids]),
        asv_count = 1
      ),
      FUN = sum
    ) |>
      transform(method = method, threshold = thresholds[i])
  }))
}))

genus_totals <- aggregate(
  flagged_reads ~ Genus,
  data = flagged_taxonomy,
  FUN = sum
)
top_genera <- head(
  genus_totals$Genus[order(genus_totals$flagged_reads, decreasing = TRUE)],
  10L
)
flagged_taxonomy$display_genus <- ifelse(
  flagged_taxonomy$Genus %in% top_genera,
  flagged_taxonomy$Genus,
  "Other"
)
taxonomy_plot_data <- aggregate(
  cbind(flagged_reads, asv_count) ~ method + threshold + display_genus,
  data = flagged_taxonomy,
  FUN = sum
)
taxonomy_plot_data$method <- factor(
  taxonomy_plot_data$method,
  levels = method_levels,
  labels = method_labels
)
taxonomy_plot_data$display_genus <- factor(
  taxonomy_plot_data$display_genus,
  levels = c(top_genera, "Other")
)
plot_taxonomy <- ggplot2::ggplot(
  taxonomy_plot_data,
  ggplot2::aes(x = factor(threshold), y = flagged_reads, fill = display_genus)
) +
  ggplot2::geom_col(width = 0.82) +
  ggplot2::facet_wrap(~method, scales = "free_y", nrow = 1) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(scale_cut = scales::cut_si(""))
  ) +
  ggplot2::scale_fill_brewer(palette = "Paired", na.value = "#cbd5e1") +
  ggplot2::labs(
    title = "Taxonomic composition of reads assigned to flagged ASVs",
    subtitle = "Top 10 genera across all runs; remaining genera grouped as Other",
    x = "decontam threshold",
    y = "Reads in flagged ASVs",
    fill = "Genus"
  ) +
  ggplot2::theme(legend.position = "right")

flag_matrix <- do.call(cbind, lapply(method_levels, function(method) {
  x <- method_results[[method]]$flags
  colnames(x) <- paste(method, colnames(x), sep = "_")
  x
}))
flagged_any <- rowSums(flag_matrix) > 0
flagged_ids <- feature_ids[flagged_any]
flagged_rank <- flagged_ids[
  order(feature_read_totals[flagged_ids], decreasing = TRUE)
]
heatmap_ids <- head(flagged_rank, 60L)
heatmap_data <- do.call(rbind, lapply(seq_along(heatmap_ids), function(i) {
  id <- heatmap_ids[i]
  labels <- colnames(flag_matrix)
  parts <- strsplit(labels, "_", fixed = TRUE)
  data.frame(
    feature_id = id,
    feature_order = i,
    method = vapply(parts, `[[`, character(1), 1L),
    threshold = as.numeric(vapply(parts, `[[`, character(1), 2L)),
    flagged = as.integer(flag_matrix[id, ]),
    stringsAsFactors = FALSE
  )
}))
heatmap_data$method <- factor(heatmap_data$method, levels = method_levels)
plot_heatmap <- ggplot2::ggplot(
  heatmap_data,
  ggplot2::aes(x = threshold, y = -feature_order, fill = factor(flagged))
) +
  ggplot2::geom_tile(width = 0.095, height = 0.94) +
  ggplot2::facet_wrap(~method, nrow = 1, labeller = ggplot2::as_labeller(method_labels)) +
  ggplot2::scale_fill_manual(
    values = c(`0` = "#e8edf2", `1` = "#155e75"),
    labels = c(`0` = "Not flagged", `1` = "Flagged")
  ) +
  ggplot2::scale_x_continuous(breaks = thresholds) +
  ggplot2::scale_y_continuous(breaks = NULL) +
  ggplot2::labs(
    title = "Flagging pattern for the 60 most abundant affected ASVs",
    subtitle = "Rows are ordered by total read abundance; each tile is one method-threshold decision",
    x = "decontam threshold",
    y = "ASVs (abundance rank)",
    fill = NULL
  )

flagged_detail <- data.frame(
  feature_id = feature_ids,
  total_reads = unname(feature_read_totals[feature_ids]),
  biological_reads = rowSums(counts[, !is_negative, drop = FALSE]),
  control_reads = rowSums(counts[, is_negative, drop = FALSE]),
  p_prevalence = prevalence$detailed[["p"]],
  p_frequency = frequency$detailed[["p"]],
  p_combined = combined$detailed[["p"]],
  prevalence_0.1 = prevalence$flags[, "0.1"],
  frequency_0.1 = frequency$flags[, "0.1"],
  combined_0.1 = combined$flags[, "0.1"],
  stringsAsFactors = FALSE
)
if (!is.null(taxonomy)) {
  flagged_detail <- cbind(
    flagged_detail,
    taxonomy[feature_ids, c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
             drop = FALSE],
    row.names = NULL
  )
}
flagged_detail <- flagged_detail[
  flagged_detail$prevalence_0.1 |
    flagged_detail$frequency_0.1 |
    flagged_detail$combined_0.1,
  ,
  drop = FALSE
]
flagged_detail <- flagged_detail[
  order(flagged_detail$total_reads, decreasing = TRUE),
  ,
  drop = FALSE
]

summary_table <- within(as.data.frame(threshold_summary), {
  method <- method_labels[as.character(method)]
  reads_removed_total_pct <- round(reads_removed_total_pct, 3)
  reads_removed_biological_pct <- round(reads_removed_biological_pct, 3)
  reads_removed_control_pct <- round(reads_removed_control_pct, 3)
})
summary_table <- summary_table[, c(
  "method", "threshold", "contaminant_features", "retained_features",
  "reads_removed_total_pct", "reads_removed_biological_pct",
  "reads_removed_control_pct"
)]
names(summary_table) <- c(
  "Method", "Threshold", "Flagged ASVs", "Retained ASVs",
  "All reads removed (%)", "Biological reads removed (%)",
  "Control reads removed (%)"
)

top_taxonomy_table <- do.call(rbind, lapply(
  split(flagged_taxonomy, interaction(flagged_taxonomy$method,
                                      flagged_taxonomy$threshold)),
  function(x) head(x[order(x$flagged_reads, decreasing = TRUE), ], 5L)
))
top_taxonomy_table$method <- method_labels[top_taxonomy_table$method]
top_taxonomy_table <- top_taxonomy_table[, c(
  "method", "threshold", "Genus", "asv_count", "flagged_reads"
)]
names(top_taxonomy_table) <- c(
  "Method", "Threshold", "Genus", "Flagged ASVs", "Flagged reads"
)
top_taxonomy_table <- top_taxonomy_table[
  order(
    match(top_taxonomy_table$Method, unname(method_labels)),
    top_taxonomy_table$Threshold,
    -top_taxonomy_table$`Flagged reads`
  ),
]

utils::write.csv(
  flagged_detail,
  file.path(output_dir, "flagged_features_threshold_0.1.csv"),
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  top_taxonomy_table,
  file.path(output_dir, "top_flagged_taxa_by_threshold.csv"),
  row.names = FALSE,
  na = ""
)

threshold_01 <- threshold_summary[threshold_summary$threshold == 0.1, ]
get_value <- function(method, column) {
  threshold_01[as.character(threshold_01$method) == method, column][[1L]]
}
concentration_cv <- 100 * stats::sd(concentration) / mean(concentration)

flagged_uri <- save_plot_uri(plot_flagged, 7.4, 4.3)
retained_uri <- save_plot_uri(plot_retained, 10.6, 4.0)
overlap_uri <- save_plot_uri(plot_overlap, 7.4, 4.3)
taxonomy_uri <- save_plot_uri(plot_taxonomy, 11.2, 5.1)
heatmap_uri <- save_plot_uri(plot_heatmap, 10.8, 5.1)

report_path <- file.path(output_dir, "decontam_threshold_comparison_v3.html")
generated_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

html <- paste0(
'<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>decontam threshold comparison - DT_Swab v3</title>
<style>
:root { --ink:#102a43; --muted:#52616f; --line:#d9e2ec; --paper:#ffffff; --wash:#f4f7fa; --accent:#007c91; --warn:#8a4b08; }
* { box-sizing:border-box; }
body { margin:0; background:var(--wash); color:var(--ink); font-family:Inter,Segoe UI,Arial,sans-serif; line-height:1.55; }
main { max-width:1180px; margin:0 auto; padding:32px 22px 64px; }
header { background:linear-gradient(135deg,#0b3543,#007c91); color:white; padding:34px; border-radius:16px; }
h1 { margin:0 0 8px; font-size:clamp(26px,4vw,42px); line-height:1.15; }
h2 { margin:38px 0 12px; font-size:24px; }
h3 { margin:20px 0 8px; font-size:18px; }
p { margin:8px 0; }
.subtitle { max-width:830px; opacity:.9; font-size:17px; }
.meta { margin-top:18px; font-size:13px; opacity:.8; }
.stats { display:grid; grid-template-columns:repeat(auto-fit,minmax(170px,1fr)); gap:12px; margin:22px 0; }
.stat { background:var(--paper); border:1px solid var(--line); border-radius:12px; padding:15px 16px; }
.stat b { display:block; font-size:27px; color:var(--accent); }
.stat span { color:var(--muted); font-size:13px; }
.note { border-left:4px solid var(--accent); background:#e7f6f8; padding:14px 16px; margin:18px 0; }
.warning { border-left-color:#d97706; background:#fff7e6; color:#633c11; }
.method-grid { display:grid; grid-template-columns:repeat(auto-fit,minmax(230px,1fr)); gap:14px; }
.method { background:var(--paper); border-top:4px solid var(--accent); padding:16px; border-radius:10px; }
.method:nth-child(2) { border-top-color:#d97706; }
.method:nth-child(3) { border-top-color:#6d5bd0; }
.figure { background:var(--paper); border:1px solid var(--line); border-radius:14px; padding:12px; margin:16px 0 24px; }
.figure img { display:block; width:100%; height:auto; }
.caption { color:var(--muted); font-size:13px; padding:4px 8px 8px; }
.table-wrap { overflow-x:auto; max-height:620px; border:1px solid var(--line); background:var(--paper); }
table { border-collapse:collapse; width:100%; font-size:13px; }
th, td { padding:8px 10px; border-bottom:1px solid var(--line); text-align:left; vertical-align:top; }
th { position:sticky; top:0; background:#eaf0f5; z-index:1; white-space:nowrap; }
tbody tr:nth-child(even) { background:#f8fafc; }
code { background:#eaf0f5; padding:2px 5px; border-radius:4px; }
footer { margin-top:38px; color:var(--muted); font-size:12px; }
@media (max-width:650px) { main{padding:14px 10px 40px;} header{padding:24px 18px;border-radius:10px;} h2{font-size:21px;} }
</style>
</head>
<body><main>
<header>
<h1>decontam threshold comparison</h1>
<div class="subtitle">Prevalence, frequency, and combined methods on the unfiltered DT_Swab <strong>v3</strong> phyloseq object. Thresholds increase from 0.1 to 0.7 in 0.1 steps.</div>
<div class="meta">Generated ', html_escape(generated_at), ' | decontam ',
html_escape(as.character(utils::packageVersion("decontam"))),
' | direct <code>decontam::isContaminant()</code> calls</div>
</header>

<section class="stats">
<div class="stat"><b>', format(length(sample_ids), big.mark = ","), '</b><span>samples</span></div>
<div class="stat"><b>', format(length(feature_ids), big.mark = ","), '</b><span>ASVs before decontam filtering</span></div>
<div class="stat"><b>', sum(is_negative), '</b><span>negative controls</span></div>
<div class="stat"><b>', format(round(concentration_cv, 1), nsmall = 1), '%</b><span>concentration CV</span></div>
</section>

<div class="note"><strong>At threshold 0.1:</strong> prevalence flagged ',
format(get_value("prevalence", "contaminant_features"), big.mark = ","),
' ASVs, frequency flagged ',
format(get_value("frequency", "contaminant_features"), big.mark = ","),
', and combined flagged ',
format(get_value("combined", "contaminant_features"), big.mark = ","),
'. Combined removed ',
format(round(get_value("combined", "reads_removed_control_pct"), 1), nsmall = 1),
'% of negative-control reads while removing ',
format(round(get_value("combined", "reads_removed_biological_pct"), 1), nsmall = 1),
'% of biological-sample reads.</div>

<div class="note warning"><strong>Frequency caution:</strong> the supplied original-QC concentrations span only ',
format(min(concentration), nsmall = 2), '-', format(max(concentration), nsmall = 2),
' ng/uL (CV ', format(round(concentration_cv, 1), nsmall = 1),
'%). The method runs correctly, but this narrow range limits frequency-based discrimination and explains the unstable jump at higher thresholds.</div>

<h2>Methods compared</h2>
<section class="method-grid">
<div class="method"><h3>Prevalence</h3><p>Uses presence/absence differences between 50 biological samples and 5 negative controls.</p></div>
<div class="method"><h3>Frequency</h3><p>Uses the inverse relationship between feature frequency and DNA concentration.</p></div>
<div class="method"><h3>Combined</h3><p>Combines prevalence and frequency evidence using decontam\'s documented combined method.</p></div>
</section>

<h2>Threshold sensitivity</h2>
<div class="figure"><img alt="Line chart of flagged ASV counts by method and threshold" src="', flagged_uri, '"><div class="caption">The decision boundary changes most sharply for prevalence at 0.6 and for frequency at 0.5-0.6.</div></div>
<div class="figure"><img alt="Faceted line chart of reads retained by sample group" src="', retained_uri, '"><div class="caption">A useful threshold should preferentially remove control reads without collapsing biological reads.</div></div>

<h2>Method agreement</h2>
<div class="figure"><img alt="Line chart of pairwise Jaccard agreement" src="', overlap_uri, '"><div class="caption">Low agreement is expected when frequency has little concentration contrast; combined is not simply the union of the two classifications.</div></div>
<div class="figure"><img alt="Heatmap of flagged ASV decisions" src="', heatmap_uri, '"><div class="caption">The heatmap exposes which abundant ASVs become flagged only after more permissive thresholds are used.</div></div>

<h2>Taxonomy of flagged reads</h2>
<div class="figure"><img alt="Stacked bar charts of flagged reads by genus" src="', taxonomy_uri, '"><div class="caption">Taxonomic composition is read-weighted; a genus can dominate because one or a few abundant ASVs are flagged.</div></div>

<h3>Top five flagged genera per run</h3>
', html_table(top_taxonomy_table, digits = 1L), '

<h2>Complete threshold summary</h2>
', html_table(summary_table, digits = 3L), '

<h2>ASVs flagged at threshold 0.1</h2>
<p>This table contains every ASV flagged by at least one method at the conventional 0.1 threshold, ordered by total reads.</p>
', html_table(flagged_detail, digits = 5L), '

<h2>Reproducibility notes</h2>
<ul>
<li>Input: <code>phy_F270R240_260812_v3.rds</code> (55 samples, 2,108 ASVs).</li>
<li>Concentration: <code>Original_QC_Conc_ng_uL</code> from the supplied sequencing metadata.</li>
<li>Negative controls: the five v3 rows with missing <code>control_status</code>, cross-checked against the Excel <code>Negative control</code> category.</li>
<li>For every method and threshold, the report verified that detailed p-values reproduce the direct logical result returned by <code>decontam::isContaminant()</code>.</li>
<li>Missing p-values are not classified as contaminants. Coverage: prevalence ',
sum(!is.na(prevalence$detailed[["p"]])), '/', length(feature_ids),
', frequency ', sum(!is.na(frequency$detailed[["p"]])), '/', length(feature_ids),
', combined ', sum(!is.na(combined$detailed[["p"]])), '/', length(feature_ids), '.</li>
</ul>
<footer>Standalone report generated from the independent example folder; it does not modify or call decontamSensitivity package code.</footer>
</main></body></html>'
)

report_connection <- file(report_path, open = "wb")
writeBin(charToRaw(enc2utf8(html)), report_connection)
close(report_connection)
cat("REPORT: ", report_path, "\n", sep = "")
