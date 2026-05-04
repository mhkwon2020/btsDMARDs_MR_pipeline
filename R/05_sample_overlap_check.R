# R/05_sample_overlap_check.R
# Sample overlap assessment utilities
# Version: 1.0
# Author: Mihye Kwon
suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

normalize_aliases <- function(x, alias_map = NULL) {
  x <- unique(trimws(as.character(x)))
  x <- x[nzchar(x)]
  if (length(x) == 0) return(character(0))

  if (!is.null(alias_map) && length(alias_map) > 0) {
    for (i in seq_along(x)) {
      key1 <- x[i]
      key2 <- tolower(x[i])
      if (!is.null(alias_map[[key1]])) {
        x[i] <- alias_map[[key1]]
      } else if (!is.null(alias_map[[key2]])) {
        x[i] <- alias_map[[key2]]
      }
    }
  }
  unique(x)
}

resolve_exposure_entry <- function(meta, exposure_id = NULL) {
  if (!is.null(meta$exposures)) {
    resolved_id <- exposure_id %||% meta$default_exposure_id
    if (is.null(resolved_id) || !resolved_id %in% names(meta$exposures)) {
      stop("Exposure id not found in metadata exposures: ", resolved_id, call. = FALSE)
    }
    entry <- meta$exposures[[resolved_id]]
    entry$exposure_id <- resolved_id
    return(entry)
  }

  if (!is.null(meta$exposure)) {
    entry <- meta$exposure
    entry$exposure_id <- exposure_id %||% "exposure"
    return(entry)
  }

  stop("Metadata must contain either `exposure:` or `exposures:`.", call. = FALSE)
}

assess_sample_overlap <- function(outcomes, metadata_file, exposure_id = NULL) {
  meta <- yaml::read_yaml(metadata_file)
  exp_entry <- resolve_exposure_entry(meta, exposure_id = exposure_id)

  exposure_aliases <- normalize_aliases(
    exp_entry$cohort_aliases,
    meta$alias_map
  )

  data.table::rbindlist(lapply(outcomes, function(o) {
    oc <- meta$outcomes[[o$name]]

    if (is.null(oc)) {
      return(data.table(
        exposure_id = exp_entry$exposure_id,
        exposure_name = exp_entry$source_name %||% exp_entry$exposure_id,
        outcome_name = o$name,
        outcome_id = o$id,
        overlap_status = "not_assessable_missing_outcome_entry",
        n_overlapping_aliases = 0L,
        overlapping_aliases = NA_character_,
        caveat_level = "high",
        recommended_action = "add outcome entry to sample_overlap_metadata.yaml"
      ))
    }

    outcome_aliases <- normalize_aliases(oc$cohort_aliases, meta$alias_map)

    if (length(outcome_aliases) == 0) {
      return(data.table(
        exposure_id = exp_entry$exposure_id,
        exposure_name = exp_entry$source_name %||% exp_entry$exposure_id,
        outcome_name = o$name,
        outcome_id = o$id,
        overlap_status = "not_assessable_missing_outcome_cohort_list",
        n_overlapping_aliases = 0L,
        overlapping_aliases = NA_character_,
        caveat_level = "high",
        recommended_action = "fill exact cohort roster from source publication"
      ))
    }

    if (length(exposure_aliases) == 0) {
      return(data.table(
        exposure_id = exp_entry$exposure_id,
        exposure_name = exp_entry$source_name %||% exp_entry$exposure_id,
        outcome_name = o$name,
        outcome_id = o$id,
        overlap_status = "not_assessable_missing_exposure_cohort_list",
        n_overlapping_aliases = 0L,
        overlapping_aliases = NA_character_,
        caveat_level = "high",
        recommended_action = "fill exact exposure cohort roster from source publication"
      ))
    }

    ov <- intersect(exposure_aliases, outcome_aliases)

    if (length(ov) == 0) {
      status <- "no_detected_direct_cohort_overlap"
      caveat <- if (!is.null(exp_entry$cohort_list_status) && identical(exp_entry$cohort_list_status, "incomplete_for_final_lock")) "moderate" else "low"
      action <- "report as no detected direct overlap based on current metadata"
    } else {
      status <- "possible_but_unquantified_overlap"
      caveat <- "moderate"
      action <- "quantify overlap if possible; otherwise retain as sensitivity caveat"
    }

    data.table(
      exposure_id = exp_entry$exposure_id,
      exposure_name = exp_entry$source_name %||% exp_entry$exposure_id,
      outcome_name = o$name,
      outcome_id = o$id,
      overlap_status = status,
      n_overlapping_aliases = length(ov),
      overlapping_aliases = if (length(ov) == 0) NA_character_ else paste(ov, collapse = "; "),
      caveat_level = caveat,
      recommended_action = action
    )
  }))
}

assess_sample_overlap_all_exposures <- function(outcomes, metadata_file, exposure_ids = NULL) {
  meta <- yaml::read_yaml(metadata_file)

  if (is.null(meta$exposures)) {
    return(assess_sample_overlap(outcomes, metadata_file, exposure_id = exposure_ids[1] %||% NULL))
  }

  ids <- exposure_ids %||% names(meta$exposures)
  data.table::rbindlist(lapply(ids, function(id) {
    assess_sample_overlap(outcomes = outcomes, metadata_file = metadata_file, exposure_id = id)
  }), fill = TRUE)
}

compute_overlap_bias_bound <- function(df, overlap_fraction = 1) {
  dt <- as.data.table(df)

  f_col <- if ("f_median" %in% names(dt)) "f_median" else if ("f_min" %in% names(dt)) "f_min" else NULL
  if (is.null(f_col)) stop("No F-statistic column found. Need f_median or f_min.")

  out <- dt[, .(
    exposure_id = if ("exposure_id" %in% names(dt)) exposure_id else NA_character_,
    gene,
    outcome_name,
    outcome_id = if ("outcome_id" %in% names(dt)) outcome_id else NA_character_,
    method = if ("method" %in% names(dt)) method else NA_character_,
    nsnp_usable = if ("nsnp_usable" %in% names(dt)) nsnp_usable else NA_integer_,
    f_statistic_used = get(f_col)
  )]

  out[, assumed_overlap_fraction := overlap_fraction]
  out[, relative_bias := assumed_overlap_fraction / f_statistic_used]
  out[, relative_bias_pct := 100 * relative_bias]
  out[, concern_level := fifelse(
    relative_bias > 0.05, "high",
    fifelse(relative_bias > 0.01, "moderate", "low")
  )]

  out[]
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}
