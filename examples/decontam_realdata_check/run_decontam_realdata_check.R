script_directory <- function() {
  command <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command, value = TRUE)
  if (length(file_arg)) {
    return(dirname(sub("^--file=", "", file_arg[1L])))
  }
  getwd()
}

analysis_dir <- script_directory()
project_dir <- file.path(analysis_dir, "..", "..")
rds_path <- file.path(
  project_dir,
  "Phyloseq",
  "phy_F270R240_260812_v3.rds"
)
concentration_path <- file.path(
  analysis_dir,
  "input",
  "concentration_metadata.csv"
)
output_dir <- file.path(analysis_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("phyloseq", "decontam")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}
if (!file.exists(rds_path)) {
  stop("Missing phyloseq RDS: ", rds_path, call. = FALSE)
}
if (!file.exists(concentration_path)) {
  stop(
    "Missing concentration CSV. Run prepare_concentration_csv.py first: ",
    concentration_path,
    call. = FALSE
  )
}

ps <- readRDS(rds_path)
if (!inherits(ps, "phyloseq")) {
  stop("The RDS does not contain a phyloseq object.", call. = FALSE)
}

metadata <- as.data.frame(phyloseq::sample_data(ps))
if (!"is.neg" %in% names(metadata)) {
  if (!"control_status" %in% names(metadata)) {
    stop(
      "phyloseq sample_data contains neither `is.neg` nor `control_status`.",
      call. = FALSE
    )
  }
  metadata[["is.neg"]] <-
    is.na(metadata[["control_status"]]) |
    metadata[["control_status"]] == "control"
}
if (anyNA(metadata[["is.neg"]])) {
  stop("`is.neg` contains missing values.", call. = FALSE)
}

concentration_data <- utils::read.csv(
  concentration_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_columns <- c(
  "Sample_Name",
  "Sample_Category",
  "Original_QC_Conc_ng_uL"
)
missing_columns <- setdiff(required_columns, names(concentration_data))
if (length(missing_columns)) {
  stop(
    "Concentration CSV is missing columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

sample_ids <- phyloseq::sample_names(ps)
excel_ids <- as.character(concentration_data[["Sample_Name"]])
if (anyDuplicated(sample_ids) || anyDuplicated(excel_ids)) {
  stop("Sample IDs must be unique in both inputs.", call. = FALSE)
}
if (!setequal(sample_ids, excel_ids)) {
  stop(
    "Sample IDs do not match exactly. RDS-only: ",
    paste(setdiff(sample_ids, excel_ids), collapse = ", "),
    "; Excel-only: ",
    paste(setdiff(excel_ids, sample_ids), collapse = ", "),
    call. = FALSE
  )
}

match_index <- match(sample_ids, excel_ids)
concentration_data <- concentration_data[match_index, , drop = FALSE]
rownames(concentration_data) <- sample_ids
concentration <- as.numeric(
  concentration_data[["Original_QC_Conc_ng_uL"]]
)
names(concentration) <- sample_ids
if (anyNA(concentration) || any(!is.finite(concentration)) ||
    any(concentration <= 0)) {
  stop(
    "Original_QC_Conc_ng_uL must be finite and positive for every sample.",
    call. = FALSE
  )
}

is_negative <- as.logical(metadata[["is.neg"]])
names(is_negative) <- rownames(metadata)
is_negative <- is_negative[sample_ids]
category_is_negative <-
  concentration_data[["Sample_Category"]] == "Negative control"
if (!identical(unname(is_negative), unname(category_is_negative))) {
  stop(
    "Excel Sample_Category does not agree with phyloseq `is.neg`.",
    call. = FALSE
  )
}

counts <- methods::as(phyloseq::otu_table(ps), "matrix")
if (!phyloseq::taxa_are_rows(ps)) counts <- t(counts)
counts <- counts[, sample_ids, drop = FALSE]
feature_ids <- rownames(counts)
thresholds <- seq(0.1, 0.7, by = 0.1)

run_method <- function(method) {
  method_args <- list(
    seqtab = ps,
    method = method,
    threshold = thresholds[1L],
    normalize = TRUE,
    detailed = TRUE
  )
  if (identical(method, "prevalence")) {
    method_args$neg <- is_negative
  } else if (identical(method, "frequency")) {
    method_args$conc <- concentration
  } else if (identical(method, "combined")) {
    method_args$neg <- is_negative
    method_args$conc <- concentration
  } else {
    stop("Unsupported method in this example: ", method, call. = FALSE)
  }

  start_time <- proc.time()[["elapsed"]]
  detailed <- do.call(decontam::isContaminant, method_args)
  flags <- vapply(thresholds, function(threshold) {
    method_args$threshold <- threshold
    method_args$detailed <- FALSE
    direct <- do.call(decontam::isContaminant, method_args)
    if (!is.logical(direct) || length(direct) != length(feature_ids) ||
        anyNA(direct)) {
      stop(
        "Direct decontam classification returned an incompatible vector for ",
        method,
        " at threshold ",
        threshold,
        ".",
        call. = FALSE
      )
    }
    if (!is.null(names(direct))) direct <- direct[feature_ids]
    unname(direct)
  }, logical(length(feature_ids)))
  rownames(flags) <- feature_ids
  colnames(flags) <- format(thresholds, trim = TRUE, scientific = FALSE)
  elapsed <- proc.time()[["elapsed"]] - start_time

  if (!is.data.frame(detailed) || !"p" %in% names(detailed) ||
      !is.numeric(detailed[["p"]]) ||
      !setequal(rownames(detailed), feature_ids)) {
    stop(
      "Detailed decontam output is incompatible for method `",
      method,
      "`.",
      call. = FALSE
    )
  }
  detailed <- detailed[feature_ids, , drop = FALSE]
  score_flags <- vapply(thresholds, function(threshold) {
    !is.na(detailed[["p"]]) & detailed[["p"]] < threshold
  }, logical(length(feature_ids)))
  if (!identical(unname(score_flags), unname(flags))) {
    stop(
      "Detailed p-values do not reproduce direct classifications for method `",
      method,
      "`.",
      call. = FALSE
    )
  }

  list(detailed = detailed, flags = flags, elapsed_seconds = elapsed)
}

prevalence <- run_method("prevalence")
frequency <- run_method("frequency")
combined <- run_method("combined")
method_results <- list(
  prevalence = prevalence,
  frequency = frequency,
  combined = combined
)

summarize_method <- function(method, result) {
  do.call(rbind, lapply(seq_along(thresholds), function(i) {
    flagged <- result$flags[, i]
    reads_total <- sum(counts)
    reads_biological <- sum(counts[, !is_negative, drop = FALSE])
    reads_control <- sum(counts[, is_negative, drop = FALSE])
    flagged_total <- sum(counts[flagged, , drop = FALSE])
    flagged_biological <- sum(
      counts[flagged, !is_negative, drop = FALSE]
    )
    flagged_control <- sum(counts[flagged, is_negative, drop = FALSE])
    data.frame(
      method = method,
      threshold = thresholds[i],
      contaminant_features = sum(flagged),
      retained_features = sum(!flagged),
      reads_removed_total = flagged_total,
      reads_removed_biological = flagged_biological,
      reads_removed_control = flagged_control,
      reads_removed_total_pct = 100 * flagged_total / reads_total,
      reads_removed_biological_pct =
        100 * flagged_biological / reads_biological,
      reads_removed_control_pct = 100 * flagged_control / reads_control,
      stringsAsFactors = FALSE
    )
  }))
}

threshold_summary <- rbind(
  summarize_method("prevalence", prevalence),
  summarize_method("frequency", frequency),
  summarize_method("combined", combined)
)

method_pairs <- utils::combn(names(method_results), 2L, simplify = FALSE)
overlap_summary <- do.call(rbind, lapply(method_pairs, function(pair) {
  do.call(rbind, lapply(seq_along(thresholds), function(i) {
    flag_a <- method_results[[pair[1L]]]$flags[, i]
    flag_b <- method_results[[pair[2L]]]$flags[, i]
    union_count <- sum(flag_a | flag_b)
    data.frame(
      threshold = thresholds[i],
      method_a = pair[1L],
      method_b = pair[2L],
      a_only = sum(flag_a & !flag_b),
      b_only = sum(!flag_a & flag_b),
      both = sum(flag_a & flag_b),
      union = union_count,
      jaccard = if (union_count == 0) NA_real_ else
        sum(flag_a & flag_b) / union_count,
      stringsAsFactors = FALSE
    )
  }))
}))

feature_results <- data.frame(
  feature_id = feature_ids,
  p_prevalence = prevalence$detailed[["p"]],
  p_frequency = frequency$detailed[["p"]],
  p_combined = combined$detailed[["p"]],
  prevalence_contaminant_0_1 = prevalence$flags[, "0.1"],
  frequency_contaminant_0_1 = frequency$flags[, "0.1"],
  combined_contaminant_0_1 = combined$flags[, "0.1"],
  stringsAsFactors = FALSE
)
taxonomy <- phyloseq::tax_table(ps, errorIfNULL = FALSE)
if (!is.null(taxonomy)) {
  taxonomy <- as.data.frame(taxonomy, stringsAsFactors = FALSE)
  taxonomy <- taxonomy[feature_ids, , drop = FALSE]
  feature_results <- cbind(feature_results, taxonomy, row.names = NULL)
}

matched_concentration <- data.frame(
  sample_id = sample_ids,
  is_negative = is_negative,
  sample_category = concentration_data[["Sample_Category"]],
  original_qc_concentration_ng_uL = concentration,
  library_size = as.numeric(phyloseq::sample_sums(ps)[sample_ids]),
  stringsAsFactors = FALSE
)

diagnostic_summary <- data.frame(
  metric = c(
    "samples_matched",
    "biological_samples",
    "negative_controls",
    "features_total",
    "prevalence_p_non_missing",
    "prevalence_p_missing",
    "frequency_p_non_missing",
    "frequency_p_missing",
    "combined_p_non_missing",
    "combined_p_missing",
    "concentration_min_ng_uL",
    "concentration_max_ng_uL",
    "concentration_mean_ng_uL",
    "concentration_cv"
  ),
  value = c(
    length(sample_ids),
    sum(!is_negative),
    sum(is_negative),
    length(feature_ids),
    sum(!is.na(prevalence$detailed[["p"]])),
    sum(is.na(prevalence$detailed[["p"]])),
    sum(!is.na(frequency$detailed[["p"]])),
    sum(is.na(frequency$detailed[["p"]])),
    sum(!is.na(combined$detailed[["p"]])),
    sum(is.na(combined$detailed[["p"]])),
    min(concentration),
    max(concentration),
    mean(concentration),
    stats::sd(concentration) / mean(concentration)
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  matched_concentration,
  file.path(output_dir, "matched_sample_concentrations.csv"),
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  threshold_summary,
  file.path(output_dir, "threshold_summary.csv"),
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  overlap_summary,
  file.path(output_dir, "method_overlap_summary.csv"),
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  feature_results,
  file.path(output_dir, "feature_level_results.csv"),
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  diagnostic_summary,
  file.path(output_dir, "diagnostic_summary.csv"),
  row.names = FALSE,
  na = ""
)
saveRDS(
  list(
    thresholds = thresholds,
    concentration = concentration,
    prevalence = prevalence,
    frequency = frequency,
    combined = combined,
    threshold_summary = threshold_summary,
    overlap_summary = overlap_summary,
    input_rds = rds_path,
    concentration_source = concentration_path
  ),
  file.path(output_dir, "decontam_realdata_check_results.rds")
)

threshold_0_1 <- threshold_summary[
  threshold_summary$threshold == 0.1,
  ,
  drop = FALSE
]
summary_lines <- c(
  "STATUS: PASS",
  paste0("Input phyloseq: ", rds_path),
  paste0("Samples: ", phyloseq::nsamples(ps)),
  paste0("Features: ", phyloseq::ntaxa(ps)),
  paste0("Biological samples: ", sum(!is_negative)),
  paste0("Negative controls: ", sum(is_negative)),
  paste0(
    "Original QC concentration range (ng/uL): ",
    format(min(concentration), trim = TRUE),
    " to ",
    format(max(concentration), trim = TRUE)
  ),
  paste0(
    "Prevalence elapsed seconds: ",
    format(round(prevalence$elapsed_seconds, 3), nsmall = 3)
  ),
  paste0(
    "Frequency elapsed seconds: ",
    format(round(frequency$elapsed_seconds, 3), nsmall = 3)
  ),
  paste0(
    "Combined elapsed seconds: ",
    format(round(combined$elapsed_seconds, 3), nsmall = 3)
  ),
  paste0("decontam version: ", as.character(utils::packageVersion("decontam"))),
  paste0(
    "Prevalence p-value coverage: ",
    sum(!is.na(prevalence$detailed[["p"]])),
    "/",
    length(feature_ids)
  ),
  paste0(
    "Frequency p-value coverage: ",
    sum(!is.na(frequency$detailed[["p"]])),
    "/",
    length(feature_ids)
  ),
  paste0(
    "Combined p-value coverage: ",
    sum(!is.na(combined$detailed[["p"]])),
    "/",
    length(feature_ids)
  ),
  paste0(
    "Concentration CV: ",
    format(
      round(100 * stats::sd(concentration) / mean(concentration), 1),
      nsmall = 1
    ),
    "%"
  ),
  "",
  "Threshold 0.1 summary:",
  paste(capture.output(print(threshold_0_1, row.names = FALSE)), collapse = "\n"),
  "",
  paste0(
    "Frequency-method caution: the original-QC concentration range is ",
    "narrow (0.40-0.56 ng/uL), so execution success does not establish ",
    "strong frequency-based discrimination."
  )
)
writeLines(summary_lines, file.path(output_dir, "run_summary.txt"))
writeLines(
  capture.output(utils::sessionInfo()),
  file.path(output_dir, "session_info.txt")
)

cat(paste(summary_lines, collapse = "\n"), "\n")
