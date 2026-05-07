# R/15_failure_mode_diagnostic.R
# Post-MR failure mode diagnostic module
# Version: 1.0
# Author: Mihye Kwon
#
# Purpose:
#   Rule-based classification of mechanistic failure modes for targets
#   labelled as Not_supported or Non_estimable in the MR benchmark.
# Configuration:
#   configs/failure_mode_diagnostic.yaml
# Controlled by:
#   run_full_study_v1.yaml -> execution_scope.run_failure_mode_diagnostic

suppressPackageStartupMessages({
  library(dplyr)
  library(httr)
  library(jsonlite)
})

# ─────────────────────────────────────────────────────────────────────────────
# R/15_failure_mode_diagnostic.R
# Post-MR Failure Mode Diagnostic Module
#
# Config:       configs/failure_mode_diagnostic.yaml
# Controlled:   run_cfg$execution_scope$run_failure_mode_diagnostic
#
#   - fetch_pipeline_internal:          reads existing pipeline CSV outputs
#   - cat7_iv_absent:                   attrition cascade + Not_supported subtyping
#   - cat3_pqtl_panel_absent:           pQTL panel presence / IV sufficiency
#   - cat5_mechanism_mismatch:          pre-specified MOA–eQTL structural mismatch
#   - cat1_mhc_ld_complexity:           MHC/xMHC locus flags
#   - cat6_post_transcriptional(partial): POST_TRANSCRIPTIONAL_MISMATCH + EQTL_PQTL_DISCORDANT
#
#   - fetch_genomic_context → cat2 (paralog/LD complexity)
#   - fetch_multiomics_qtl  → cat4 (multi-omics) + cat6 full (sQTL/mQTL)
#   - fetch_pharmacology    → cat5 API validation
# ─────────────────────────────────────────────────────────────────────────────


# ─────────────────────────────────────────────────────────────────────────────
# INTERNAL HELPERS
# ─────────────────────────────────────────────────────────────────────────────

failmode_read_csv <- function(path) {
  if (!fs::file_exists(path)) return(data.frame())
  tryCatch(
    read.csv(as.character(path), stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) {
      emsg <- conditionMessage(e)
      warning(glue("[fail_mode] Cannot read {path}: {emsg}"))
      data.frame()
    }
  )
}

failmode_get_int <- function(df_row, col) {
  if (is.null(df_row) || nrow(df_row) == 0 || !col %in% names(df_row)) return(NA_integer_)
  suppressWarnings(as.integer(df_row[[col]][1]))
}

failmode_get_num <- function(df_row, col) {
  if (is.null(df_row) || nrow(df_row) == 0 || !col %in% names(df_row)) return(NA_real_)
  suppressWarnings(as.numeric(df_row[[col]][1]))
}

failmode_get_chr <- function(df_row, col) {
  if (is.null(df_row) || nrow(df_row) == 0 || !col %in% names(df_row)) return(NA_character_)
  as.character(df_row[[col]][1])
}

write_csv_safe <- function(df, path) {
  tryCatch({
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(df, file = path, row.names = FALSE)
  }, error = function(e) warning("[write_csv_safe] ", conditionMessage(e)))
  invisible(path)
}

failmode_log <- function(log_file, msg, ...) {
  if (!is.null(log_file) && is.character(log_file)) {
    cat(sprintf(msg, ...), file = log_file, append = TRUE)
  }
}


# ─────────────────────────────────────────────────────────────────────────────
# CONFIG LOADER
# ─────────────────────────────────────────────────────────────────────────────

load_fail_mode_config <- function(config_path) {
  if (!fs::file_exists(config_path)) {
    stop(glue("Failure mode diagnostic config not found: {config_path}"))
  }
  yaml::read_yaml(config_path)
}


# ─────────────────────────────────────────────────────────────────────────────
# FETCH LAYER: pipeline_internal
#
# ─────────────────────────────────────────────────────────────────────────────

fetch_pipeline_internal <- function(results_dir, log_file = NULL) {
  failmode_log(log_file,
    "[fail_mode] fetch_pipeline_internal start results_dir={results_dir}")

  cache <- list(
    classification       = failmode_read_csv(
                             fs::path(results_dir, "benchmark_classification_table.csv")),
    classification_sup   = failmode_read_csv(
                             fs::path(results_dir, "benchmark_classification_with_support_layers.csv")),
    targets              = failmode_read_csv(
                             fs::path(results_dir, "master_benchmark_table.csv")),
    # Attrition cascade (S0–S4)
    attrition            = failmode_read_csv(
                             fs::path(results_dir, "instrument_attrition_table.csv")),
    mr_results           = failmode_read_csv(
                             fs::path(results_dir, "preferred_mr_results.csv")),
    # Colocalisation
    coloc                = failmode_read_csv(
                             fs::path(results_dir, "priority_target_colocalisation_summary.csv")),
    # pQTL
    pqtl_results         = failmode_read_csv(
                             fs::path(results_dir, "pqtl_main_mr_results.csv")),
    pqtl_preferred       = failmode_read_csv(
                             fs::path(results_dir, "pqtl_preferred_mr_results.csv")),
    pqtl_estimability    = failmode_read_csv(
                             fs::path(results_dir, "pqtl_non_estimable_table.csv")),
    pqtl_source_coverage = failmode_read_csv(
                             fs::path(results_dir, "pqtl_source_summary.csv")),
    pqtl_support         = failmode_read_csv(
                             fs::path(results_dir, "pqtl_support_summary.csv")),
    lead_variants        = failmode_read_csv(
                             fs::path(results_dir, "benchmark_with_lead_variant.csv"))
  )

  n_loaded <- sum(vapply(cache, function(df) !is.null(df) && nrow(df) > 0, logical(1)))
  failmode_log(log_file,
    "[fail_mode] fetch_pipeline_internal done tables_with_data={n_loaded}/{length(cache)}")

  cache
}


# ─────────────────────────────────────────────────────────────────────────────
# DIAGNOSTIC TARGET IDENTIFICATION
#
# ─────────────────────────────────────────────────────────────────────────────

identify_diagnostic_targets <- function(internal_cache, diag_cfg) {
  scope    <- diag_cfg$diagnostic_scope %||% list()
  allowed  <- scope$target_classes %||% c("Not_supported", "Non_estimable")
  if (isTRUE(scope$include_weak_or_partial)) allowed <- union(allowed, "Weak_or_partial")
  if (isTRUE(scope$exclude_strong))          allowed <- setdiff(allowed, "Strong")

  df <- internal_cache$classification
  if (is.null(df) || nrow(df) == 0 || !"class" %in% names(df)) {
    df <- internal_cache$targets
  }
  if (is.null(df) || nrow(df) == 0 || !"gene" %in% names(df)) return(character(0))

  if ("class" %in% names(df)) {
    df <- df[df$class %in% allowed, , drop = FALSE]
  }

  sort(unique(tolower(df$gene)))
}


# ─────────────────────────────────────────────────────────────────────────────
#
# Returns: list(flag, confidence, attrition_stage)
# ─────────────────────────────────────────────────────────────────────────────

classify_iv_attrition <- function(attr_row) {
  empty <- list(flag = NA_character_, confidence = NA_character_, attrition_stage = NA_integer_)
  if (is.null(attr_row) || nrow(attr_row) == 0) return(empty)

  s0 <- failmode_get_int(attr_row, "raw_variant_count")
  s1 <- failmode_get_int(attr_row, "post_p_threshold_count")
  s2 <- failmode_get_int(attr_row, "post_clump_count")
  s3 <- failmode_get_int(attr_row, "post_outcome_extraction_count")
  s4 <- failmode_get_int(attr_row, "post_harmonisation_usable_count")

  if (!is.na(s0) && s0 == 0L) {
    return(list(flag = "IV_NO_EQTL_DATA",
                confidence = "High", attrition_stage = 0L))
  }
  if (!is.na(s0) && s0 > 0L && !is.na(s1) && s1 == 0L) {
    return(list(flag = "IV_NO_SIGNIFICANT_VARIANT",
                confidence = "High", attrition_stage = 1L))
  }
  if (!is.na(s1) && s1 > 0L && !is.na(s2) && s2 == 0L) {
    return(list(flag = "IV_CLUMPING_LOSS",
                confidence = "Moderate", attrition_stage = 2L))
  }
  if (!is.na(s2) && s2 > 0L && !is.na(s3) && s3 == 0L) {
    return(list(flag = "IV_OUTCOME_MISMATCH",
                confidence = "Moderate", attrition_stage = 3L))
  }
  if (!is.na(s3) && s3 > 0L && !is.na(s4) && s4 == 0L) {
    return(list(flag = "IV_HARMONISATION_FAILURE",
                confidence = "Moderate", attrition_stage = 4L))
  }

  list(flag = NA_character_, confidence = NA_character_, attrition_stage = 5L)
}


# ─────────────────────────────────────────────────────────────────────────────
# CAT7 HELPER: Not_supported / Weak_or_partial subtype classification
#
# Returns: list(subtype, confidence)
# ─────────────────────────────────────────────────────────────────────────────

classify_not_supported_subtype_cat7 <- function(gene, class_df, coloc_df) {
  empty <- list(subtype = NA_character_, confidence = NA_character_)

  if (is.null(class_df) || nrow(class_df) == 0) return(empty)
  row <- class_df[tolower(class_df$gene) == tolower(gene), , drop = FALSE]
  if (nrow(row) == 0) return(empty)

  cls         <- failmode_get_chr(row, "class")
  support_sub <- failmode_get_chr(row, "support_subtype")
  primary_nom <- if ("primary_nominal_support" %in% names(row)) row$primary_nominal_support[1] else NA
  rep_nom     <- if ("replication_nominal_support" %in% names(row)) row$replication_nominal_support[1] else NA

  if (!is.na(cls) && cls == "Non_estimable") return(empty)

  if (!is.na(cls) && cls == "Not_supported") {
    if (!is.na(support_sub) && support_sub == "inverse") {
      return(list(subtype = "NS_INVERSE_DIRECTION", confidence = "High"))
    }

    coloc_h4 <- NA_real_
    if (!is.null(coloc_df) && nrow(coloc_df) > 0 && "gene" %in% names(coloc_df)) {
      crow <- coloc_df[tolower(coloc_df$gene) == tolower(gene), , drop = FALSE]
      if (nrow(crow) > 0) {
        coloc_h4 <- if ("coloc_pp_h4" %in% names(crow)) {
          failmode_get_num(crow, "coloc_pp_h4")
        } else if ("final_PP.H4" %in% names(crow)) {
          failmode_get_num(crow, "final_PP.H4")
        } else {
          NA_real_
        }
      }
    }

    if (isTRUE(primary_nom) && !is.na(coloc_h4) && is.finite(coloc_h4) && coloc_h4 < 0.8) {
      return(list(subtype = "NS_COLOC_NON_SUPPORT", confidence = "High"))
    }

    if (isTRUE(primary_nom) && identical(rep_nom, FALSE)) {
      return(list(subtype = "NS_NO_REPLICATION", confidence = "Moderate"))
    }

    return(list(subtype = "NS_NULL_PRIMARY", confidence = "High"))
  }

  if (!is.na(cls) && cls == "Weak_or_partial") {
    return(list(subtype = "WP_PARTIAL_SUPPORT", confidence = "Moderate"))
  }

  empty
}


# ─────────────────────────────────────────────────────────────────────────────
# CAT7 HELPER: Brion et al. (2013) power estimate
#
#
# ─────────────────────────────────────────────────────────────────────────────

estimate_mr_power_cat7 <- function(mr_se, pow_cfg) {
  if (!isTRUE(pow_cfg$enabled %||% FALSE)) {
    return(list(power_estimate = NA_real_, power_sufficient = NA))
  }
  if (is.na(mr_se) || !is.finite(mr_se) || mr_se <= 0) {
    return(list(power_estimate = NA_real_, power_sufficient = NA))
  }

  assumed_or <- pow_cfg$assumed_or  %||% 1.1
  alpha      <- pow_cfg$alpha       %||% 0.05
  thresh     <- pow_cfg$power_sufficient_threshold %||% 0.80

  z_alpha   <- qnorm(1 - alpha / 2)
  true_b    <- abs(log(assumed_or))
  power_est <- pnorm(true_b / mr_se - z_alpha)

  list(
    power_estimate   = round(power_est, 4),
    power_sufficient = power_est >= thresh
  )
}


# ─────────────────────────────────────────────────────────────────────────────
# ALTERNATIVE STRATEGY LOOKUP
# ─────────────────────────────────────────────────────────────────────────────

lookup_alternative_strategy <- function(flag, diag_cfg) {
  if (is.na(flag) || is.null(flag) || !nzchar(flag)) return(NA_character_)
  rec <- diag_cfg$alternative_strategy$recommendations[[flag]]
  if (is.null(rec)) return(NA_character_)
  rec$strategy %||% NA_character_
}


# ─────────────────────────────────────────────────────────────────────────────
# CATEGORY 7: IV Absent — main function
#
#   - Non_estimable: attrition cascade → IV absence flag
# ─────────────────────────────────────────────────────────────────────────────

run_cat7_iv_absent <- function(internal_cache, diag_cfg, log_file = NULL) {
  cfg7 <- diag_cfg$cat7_iv_absent %||% list()
  if (!isTRUE(cfg7$enabled %||% TRUE)) {
    failmode_log(log_file, "[fail_mode] cat7_iv_absent disabled in config, skipping")
    return(data.frame())
  }

  target_genes <- identify_diagnostic_targets(internal_cache, diag_cfg)
  if (length(target_genes) == 0) {
    failmode_log(log_file, "[fail_mode] cat7: no diagnostic target genes found")
    return(data.frame())
  }

  failmode_log(log_file,
    "[fail_mode] cat7 start n={length(target_genes)} genes={paste(target_genes, collapse=';')}")

  attrition_df <- internal_cache$attrition
  class_df     <- internal_cache$classification
  mr_df        <- internal_cache$mr_results
  coloc_df     <- internal_cache$coloc
  pow_cfg      <- cfg7$power_estimate %||% list()

  rows <- lapply(target_genes, function(g) {

    # ── 1. classification row ──────────────────────────────────────────────
    class_row <- if (!is.null(class_df) && "gene" %in% names(class_df)) {
      class_df[tolower(class_df$gene) == g, , drop = FALSE]
    } else {
      data.frame()
    }
    gene_class <- failmode_get_chr(class_row, "class")

    attr_all <- if (!is.null(attrition_df) && "gene" %in% names(attrition_df)) {
      attrition_df[tolower(attrition_df$gene) == g, , drop = FALSE]
    } else {
      data.frame()
    }

    attr_primary <- if (nrow(attr_all) > 0 && "outcome_name" %in% names(attr_all)) {
      r <- attr_all[attr_all$outcome_name == "RA_EUR_PRIMARY", , drop = FALSE]
      if (nrow(r) > 0) r[1, , drop = FALSE] else attr_all[1, , drop = FALSE]
    } else if (nrow(attr_all) > 0) {
      attr_all[1, , drop = FALSE]
    } else {
      data.frame()
    }

    non_est_reason <- failmode_get_chr(attr_primary, "non_estimable_reason")

    s0 <- failmode_get_int(attr_primary, "raw_variant_count")
    s1 <- failmode_get_int(attr_primary, "post_p_threshold_count")
    s2 <- failmode_get_int(attr_primary, "post_clump_count")
    s3 <- failmode_get_int(attr_primary, "post_outcome_extraction_count")
    s4 <- failmode_get_int(attr_primary, "post_harmonisation_usable_count")

    is_non_estimable <- isTRUE(!is.na(gene_class) && gene_class == "Non_estimable") ||
      (nrow(attr_primary) > 0 &&
       "estimable" %in% names(attr_primary) &&
       identical(attr_primary$estimable[1], FALSE))

    iv_result <- if (is_non_estimable) {
      classify_iv_attrition(attr_primary)
    } else {
      list(flag = NA_character_, confidence = NA_character_, attrition_stage = 5L)
    }

    ns_result <- if (!is_non_estimable) {
      classify_not_supported_subtype_cat7(g, class_df, coloc_df)
    } else {
      list(subtype = NA_character_, confidence = NA_character_)
    }

    cat7_flag <- if (!is.na(iv_result$flag)) {
      iv_result$flag
    } else if (!is.na(ns_result$subtype)) {
      ns_result$subtype
    } else {
      NA_character_
    }

    cat7_confidence <- if (!is.na(iv_result$flag)) {
      iv_result$confidence
    } else if (!is.na(ns_result$subtype)) {
      ns_result$confidence
    } else {
      NA_character_
    }

    # ── 7. Preferred MR result (RA_EUR_PRIMARY) ───────────────────────────
    mr_primary <- if (!is.null(mr_df) && nrow(mr_df) > 0 && "gene" %in% names(mr_df)) {
      r <- mr_df[tolower(mr_df$gene) == g, , drop = FALSE]
      if ("outcome_name" %in% names(r)) {
        r2 <- r[r$outcome_name == "RA_EUR_PRIMARY", , drop = FALSE]
        if (nrow(r2) > 0) r2[1, , drop = FALSE] else if (nrow(r) > 0) r[1, , drop = FALSE] else data.frame()
      } else if (nrow(r) > 0) {
        r[1, , drop = FALSE]
      } else {
        data.frame()
      }
    } else {
      data.frame()
    }

    mr_se        <- failmode_get_num(mr_primary, "se")
    mr_pval      <- failmode_get_num(mr_primary, "pval")
    mr_beta      <- failmode_get_num(mr_primary, "b")
    f_median     <- failmode_get_num(mr_primary, "f_median")
    n_f_gt_10    <- failmode_get_int(mr_primary, "n_f_gt_10")
    nsnp_usable  <- failmode_get_int(mr_primary, "nsnp_usable")

    # ── 8. Power estimate ─────────────────────────────────────────────────
    pow_result <- estimate_mr_power_cat7(mr_se, pow_cfg)

    # ── 9. Alternative strategy ───────────────────────────────────────────
    alt_strategy <- lookup_alternative_strategy(cat7_flag, diag_cfg)

    data.frame(
      gene                         = g,
      class                        = gene_class,
      cat7_flag                    = cat7_flag,
      cat7_confidence              = cat7_confidence,
      cat7_is_non_estimable        = is_non_estimable,
      cat7_ns_subtype              = ns_result$subtype,
      cat7_attrition_stage         = iv_result$attrition_stage,
      non_estimable_reason         = non_est_reason,
      attrition_s0_raw             = s0,
      attrition_s1_post_p          = s1,
      attrition_s2_post_clump      = s2,
      attrition_s3_post_extraction = s3,
      attrition_s4_usable          = s4,
      mr_beta                      = mr_beta,
      mr_se                        = mr_se,
      mr_pval                      = mr_pval,
      f_median                     = f_median,
      n_f_gt_10                    = n_f_gt_10,
      nsnp_usable                  = nsnp_usable,
      power_estimate               = pow_result$power_estimate,
      power_sufficient             = pow_result$power_sufficient,
      recommended_strategy         = alt_strategy,
      annotation_timestamp         = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
      stringsAsFactors             = FALSE
    )
  })

  result <- dplyr::bind_rows(rows)

  n_flagged <- sum(!is.na(result$cat7_flag), na.rm = TRUE)
  failmode_log(log_file,
    "[fail_mode] cat7 done n_genes={nrow(result)} n_flagged={n_flagged}")

  result
}


# ─────────────────────────────────────────────────────────────────────────────
# TIER A HELPERS
# ─────────────────────────────────────────────────────────────────────────────

get_target_annotation <- function(gene, field, diag_cfg) {
  annots <- diag_cfg$target_annotations %||% list()
  val    <- annots[[tolower(gene)]][[field]]
  if (is.null(val)) NA_character_ else as.character(val[1])
}

get_pqtl_panel_source_ids <- function(diag_cfg) {
  cfg_path <- diag_cfg$fetch_pipeline_internal$pqtl_panel_manifest$config_file %||%
    "configs/pqtl_sources.yaml"
  if (!fs::file_exists(cfg_path)) return(character(0))
  pqtl_cfg <- tryCatch(yaml::read_yaml(cfg_path), error = function(e) list())
  sources  <- pqtl_cfg$sources %||% list()
  ids      <- vapply(sources, function(s) s$source_id %||% NA_character_, character(1))
  enabled  <- vapply(sources, function(s) isTRUE(s$enabled %||% TRUE), logical(1))
  ids[enabled & !is.na(ids)]
}


# ─────────────────────────────────────────────────────────────────────────────
#
# ─────────────────────────────────────────────────────────────────────────────

run_cat3_pqtl_panel_absent <- function(internal_cache, diag_cfg, log_file = NULL) {
  cfg3 <- diag_cfg$cat3_pqtl_panel_absent %||% list()
  if (!isTRUE(cfg3$enabled %||% TRUE)) {
    failmode_log(log_file, "[fail_mode] cat3_pqtl_panel_absent disabled, skipping")
    return(data.frame())
  }

  target_genes <- identify_diagnostic_targets(internal_cache, diag_cfg)
  if (length(target_genes) == 0) return(data.frame())

  failmode_log(log_file, "[fail_mode] cat3 start n={length(target_genes)}")

  pqtl_pref    <- internal_cache$pqtl_preferred
  pqtl_nonest  <- internal_cache$pqtl_estimability
  pqtl_support <- internal_cache$pqtl_support
  class_df     <- internal_cache$classification

  rows <- lapply(target_genes, function(g) {
    # class
    class_row  <- if (!is.null(class_df) && "gene" %in% names(class_df))
      class_df[tolower(class_df$gene) == g, , drop = FALSE] else data.frame()
    gene_class <- failmode_get_chr(class_row, "class")

    pref_rows <- if (!is.null(pqtl_pref) && "gene" %in% names(pqtl_pref))
      pqtl_pref[tolower(pqtl_pref$gene) == g, , drop = FALSE] else data.frame()

    nonest_rows <- if (!is.null(pqtl_nonest) && "gene" %in% names(pqtl_nonest))
      pqtl_nonest[tolower(pqtl_nonest$gene) == g, , drop = FALSE] else data.frame()

    sup_row <- if (!is.null(pqtl_support) && "gene" %in% names(pqtl_support))
      pqtl_support[tolower(pqtl_support$gene) == g, , drop = FALSE] else data.frame()

    pqtl_estimable_any <- if (nrow(sup_row) > 0 && "pqtl_estimable_any" %in% names(sup_row))
      isTRUE(sup_row$pqtl_estimable_any[1])
    else
      nrow(pref_rows) > 0

    pqtl_n_sources <- if (nrow(sup_row) > 0 && "pqtl_n_sources" %in% names(sup_row))
      suppressWarnings(as.integer(sup_row$pqtl_n_sources[1]))
    else
      NA_integer_

    # target_annotations fallback
    annot_panel_status <- get_target_annotation(g, "pqtl_panel_status", diag_cfg)

    result <- if (pqtl_estimable_any) {
      list(flag = "PQTL_AVAILABLE", confidence = NA_character_)

    } else if (nrow(nonest_rows) > 0) {
      rv <- suppressWarnings(as.integer(nonest_rows$raw_variant_count[1]))
      reason <- if ("non_estimable_reason" %in% names(nonest_rows))
        nonest_rows$non_estimable_reason[1] else NA_character_
      source_unavailable <- !is.na(rv) && rv == 0L &&
        !is.na(reason) && reason == "no_variant_after_target_filter"
      if (source_unavailable)
        list(flag = "PQTL_PANEL_ABSENT",    confidence = "High")
      else
        list(flag = "PQTL_IV_INSUFFICIENT", confidence = "Moderate")

    } else {
      if (!is.na(annot_panel_status) && annot_panel_status == "available")
        list(flag = "PQTL_IV_INSUFFICIENT", confidence = "Low")
      else if (!is.na(annot_panel_status) && annot_panel_status == "check_required")
        list(flag = "PQTL_PANEL_ABSENT",    confidence = "Low")
      else
        list(flag = "PQTL_PANEL_ABSENT",    confidence = "Low")
    }

    alt_strategy <- lookup_alternative_strategy(result$flag, diag_cfg)

    data.frame(
      gene                  = g,
      class                 = gene_class,
      cat3_flag             = result$flag,
      cat3_confidence       = result$confidence,
      pqtl_estimable_any    = pqtl_estimable_any,
      pqtl_n_sources        = pqtl_n_sources,
      pqtl_panel_status_annot = annot_panel_status,
      recommended_strategy  = alt_strategy,
      annotation_timestamp  = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
      stringsAsFactors      = FALSE
    )
  })

  result <- dplyr::bind_rows(rows)
  n_flagged <- sum(!is.na(result$cat3_flag) & result$cat3_flag != "PQTL_AVAILABLE",
                   na.rm = TRUE)
  failmode_log(log_file,
    "[fail_mode] cat3 done n_genes={nrow(result)} n_flagged_absent_or_insufficient={n_flagged}")
  result
}


# ─────────────────────────────────────────────────────────────────────────────
#
# ─────────────────────────────────────────────────────────────────────────────

run_cat5_mechanism_mismatch <- function(internal_cache, diag_cfg, log_file = NULL,
                                        pharmacology_cache = NULL) {
  cfg5 <- diag_cfg$cat5_mechanism_mismatch %||% list()
  if (!isTRUE(cfg5$enabled %||% TRUE)) {
    failmode_log(log_file, "[fail_mode] cat5_mechanism_mismatch disabled, skipping")
    return(data.frame())
  }

  target_genes <- identify_diagnostic_targets(internal_cache, diag_cfg)
  if (length(target_genes) == 0) return(data.frame())

  failmode_log(log_file, "[fail_mode] cat5 start n={length(target_genes)}")

  class_df <- internal_cache$classification

  mismatch_rules <- cfg5$mismatch_rules %||% list()
  enzymatic_kw  <- tolower(unlist(mismatch_rules$enzymatic_inhibition$drug_moa_keywords %||%
    c("kinase inhibitor", "enzyme inhibitor", "ATP-competitive", "allosteric inhibitor")))
  secreted_kw   <- tolower(unlist(mismatch_rules$secreted_neutralisation$drug_moa_keywords %||%
    c("neutralising antibody", "ligand sequestration", "anti-cytokine", "TNF inhibitor")))
  depletion_kw  <- tolower(unlist(mismatch_rules$cell_depletion$drug_moa_keywords %||%
    c("cell depletion", "ADCC", "CDC", "anti-CD20")))
  costim_kw     <- tolower(unlist(mismatch_rules$costimulation_blockade$drug_moa_keywords %||%
    c("co-stimulation blocker", "fusion protein", "CTLA4-Ig")))
  receptor_kw   <- tolower(unlist(mismatch_rules$receptor_blockade_complex$drug_moa_keywords %||%
    c("receptor blockade", "receptor antagonism", "anti-receptor antibody")))

  keyword_match_flag <- function(ot_mech_str) {
    if (is.na(ot_mech_str)) return(NA_character_)
    m <- tolower(ot_mech_str)
    if (any(vapply(enzymatic_kw, function(k) grepl(k, m, fixed = TRUE), logical(1))))
      return("MECHANISM_MISMATCH_ENZYMATIC")
    if (any(vapply(secreted_kw,  function(k) grepl(k, m, fixed = TRUE), logical(1))))
      return("MECHANISM_MISMATCH_SECRETED")
    if (any(vapply(depletion_kw, function(k) grepl(k, m, fixed = TRUE), logical(1))))
      return("MECHANISM_MISMATCH_DEPLETION")
    if (any(vapply(costim_kw,    function(k) grepl(k, m, fixed = TRUE), logical(1))))
      return("MECHANISM_MISMATCH_COSTIM")
    if (any(vapply(receptor_kw,  function(k) grepl(k, m, fixed = TRUE), logical(1))))
      return("MECHANISM_MISMATCH_RECEPTOR_COMPLEX")
    NA_character_
  }

  rows <- lapply(target_genes, function(g) {
    class_row  <- if (!is.null(class_df) && "gene" %in% names(class_df))
      class_df[tolower(class_df$gene) == g, , drop = FALSE] else data.frame()
    gene_class <- failmode_get_chr(class_row, "class")

    # ── 1. Pre-specified annotation (primary source) ───────────────────────
    expected_mismatch <- get_target_annotation(g, "expected_mismatch",  diag_cfg)
    drug_modality     <- get_target_annotation(g, "drug_modality",      diag_cfg)
    drug_action_type  <- get_target_annotation(g, "drug_action_type",   diag_cfg)
    subcellular_loc   <- get_target_annotation(g, "subcellular_location", diag_cfg)

    pre_flag <- if (!is.na(expected_mismatch) && nzchar(expected_mismatch))
      expected_mismatch else NA_character_

    api_flag           <- NA_character_
    api_subcell        <- NA_character_
    api_source         <- NA_character_
    ot_mechanisms_str  <- NA_character_

    if (!is.null(pharmacology_cache) && !is.null(pharmacology_cache[[g]])) {
      pc               <- pharmacology_cache[[g]]
      ot_mechanisms_str <- pc$ot_mechanisms %||% NA_character_
      api_subcell      <- pc$uniprot_subcellular_location %||% NA_character_

      # Keyword match on OT action types
      api_flag_kw <- keyword_match_flag(ot_mechanisms_str)

      # If subcellular location is "Secreted" and no stronger flag from keywords, infer
      if (is.na(api_flag_kw) && !is.na(api_subcell) &&
          grepl("secreted", api_subcell, ignore.case = TRUE)) {
        api_flag_kw <- "MECHANISM_MISMATCH_SECRETED"
      }
      api_flag   <- api_flag_kw
      api_source <- if (!is.na(api_flag)) "open_targets_api" else NA_character_
    }

    # ── 3. Final flag: pre-specified is primary; OT cross-validates ──────────
    # Even when pre_flag is set, run OT keyword match and compare.
    # Discordance (both non-NA but different) downgrades confidence to "Moderate"
    # and records a note — it does NOT override the pre-specified flag.
    cat5_flag <- if (!is.na(pre_flag)) pre_flag else api_flag

    cross_validation       <- NA_character_   # "concordant" | "discordant" | "api_only" | "pre_only" | "no_api_data"
    cross_validation_note  <- NA_character_
    cross_validation_concordant <- NA          # logical

    if (!is.na(pre_flag) && !is.na(api_flag)) {
      if (pre_flag == api_flag) {
        cross_validation             <- "concordant"
        cross_validation_concordant  <- TRUE
        cross_validation_note        <- paste0("OT confirms: ", api_flag)
      } else {
        cross_validation             <- "discordant"
        cross_validation_concordant  <- FALSE
        cross_validation_note        <- paste0("pre=", pre_flag, " vs OT=", api_flag)
      }
    } else if (!is.na(pre_flag) && is.na(api_flag)) {
      cross_validation             <- if (!is.null(pharmacology_cache) &&
                                           !is.null(pharmacology_cache[[g]])) "pre_only" else "no_api_data"
      cross_validation_concordant  <- NA
    } else if (is.na(pre_flag) && !is.na(api_flag)) {
      cross_validation             <- "api_only"
      cross_validation_concordant  <- NA
    }

    # Confidence: High when concordant or pre-specified-only (no OT data); Moderate otherwise
    cat5_confidence <- if (!is.na(pre_flag)) {
      if (identical(cross_validation, "discordant")) "Moderate" else "High"
    } else if (!is.na(api_flag)) "Moderate" else NA_character_

    annotation_source <- if (!is.na(pre_flag)) "pre_specified"
    else if (!is.na(api_flag)) api_source
    else NA_character_

    alt_strategy <- lookup_alternative_strategy(cat5_flag, diag_cfg)

    data.frame(
      gene                        = g,
      class                       = gene_class,
      cat5_flag                   = cat5_flag,
      cat5_confidence             = cat5_confidence,
      drug_modality               = drug_modality,
      drug_action_type            = drug_action_type,
      subcellular_location        = if (!is.na(api_subcell)) api_subcell else subcellular_loc,
      ot_mechanisms               = ot_mechanisms_str,
      annotation_source           = annotation_source,
      cat5_cross_validation       = cross_validation,
      cat5_cross_validation_note  = cross_validation_note,
      cross_validation_concordant = cross_validation_concordant,
      recommended_strategy        = alt_strategy,
      annotation_timestamp        = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
      stringsAsFactors            = FALSE
    )
  })

  result <- dplyr::bind_rows(rows)
  n_flagged <- sum(!is.na(result$cat5_flag), na.rm = TRUE)
  failmode_log(log_file,
    "[fail_mode] cat5 done n_genes={nrow(result)} n_flagged={n_flagged}")
  result
}


# ─────────────────────────────────────────────────────────────────────────────
#
# ─────────────────────────────────────────────────────────────────────────────

run_cat1_mhc_ld_complexity <- function(internal_cache, diag_cfg, log_file = NULL) {
  cfg1 <- diag_cfg$cat1_mhc_ld_complexity %||% list()
  if (!isTRUE(cfg1$enabled %||% TRUE)) {
    failmode_log(log_file, "[fail_mode] cat1_mhc_ld_complexity disabled, skipping")
    return(data.frame())
  }

  target_genes <- identify_diagnostic_targets(internal_cache, diag_cfg)
  if (length(target_genes) == 0) return(data.frame())

  failmode_log(log_file, "[fail_mode] cat1 start n={length(target_genes)}")

  mhc_flags <- c("MHC_CORE", "MHC_EXTENDED")

  known_clusters <- diag_cfg$fetch_genomic_context$paralog_query$known_clusters %||% list()
  mhc_cluster_genes <- character(0)
  for (cluster in known_clusters) {
    if (!is.null(cluster$chromosome) && cluster$chromosome == "6") {
      mhc_cluster_genes <- union(mhc_cluster_genes, tolower(cluster$genes %||% character(0)))
    }
  }

  class_df <- internal_cache$classification

  rows <- lapply(target_genes, function(g) {
    class_row  <- if (!is.null(class_df) && "gene" %in% names(class_df))
      class_df[tolower(class_df$gene) == g, , drop = FALSE] else data.frame()
    gene_class <- failmode_get_chr(class_row, "class")

    locus_flags_raw <- diag_cfg$target_annotations[[g]]$known_locus_flags
    locus_flags     <- if (!is.null(locus_flags_raw)) tolower(locus_flags_raw) else character(0)

    mhc_hits <- locus_flags[locus_flags %in% tolower(mhc_flags)]

    in_mhc_cluster <- g %in% mhc_cluster_genes

    cat1_flag <- if (length(mhc_hits) > 0) {
      if ("mhc_core" %in% mhc_hits) "MHC_CORE" else "MHC_EXTENDED"
    } else if (in_mhc_cluster) {
      "MHC_EXTENDED"
    } else {
      NA_character_
    }

    cat1_confidence <- if (!is.na(cat1_flag)) "High" else NA_character_
    annotation_source <- if (!is.na(cat1_flag)) "pre_specified" else NA_character_

    alt_strategy <- lookup_alternative_strategy(cat1_flag, diag_cfg)

    data.frame(
      gene                 = g,
      class                = gene_class,
      cat1_flag            = cat1_flag,
      cat1_confidence      = cat1_confidence,
      known_locus_flags    = if (length(locus_flags) > 0) paste(locus_flags, collapse = "|") else NA_character_,
      in_mhc_cluster       = in_mhc_cluster,
      annotation_source    = annotation_source,
      recommended_strategy = alt_strategy,
      annotation_timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
      stringsAsFactors     = FALSE
    )
  })

  result <- dplyr::bind_rows(rows)
  n_flagged <- sum(!is.na(result$cat1_flag), na.rm = TRUE)
  failmode_log(log_file,
    "[fail_mode] cat1 done n_genes={nrow(result)} n_flagged={n_flagged}")
  result
}


# ─────────────────────────────────────────────────────────────────────────────
#
#   POST_TRANSCRIPTIONAL_MISMATCH — target_annotations.post_transcriptional_note
#
#   SPLICING_QTL_DETECTED    — eQTL Catalogue sQTL
#   ALTERNATIVE_QTL_AVAILABLE — GoDMC mQTL / QTLbase2
# ─────────────────────────────────────────────────────────────────────────────

run_cat6_post_transcriptional <- function(internal_cache, diag_cfg, log_file = NULL,
                                           qtl_cache = NULL) {
  cfg6 <- diag_cfg$cat6_post_transcriptional %||% list()
  if (!isTRUE(cfg6$enabled %||% TRUE)) {
    failmode_log(log_file, "[fail_mode] cat6_post_transcriptional disabled, skipping")
    return(data.frame())
  }

  target_genes <- identify_diagnostic_targets(internal_cache, diag_cfg)
  if (length(target_genes) == 0) return(data.frame())

  tier_label <- if (!is.null(qtl_cache)) "full" else "partial"
  failmode_log(log_file,
    "[fail_mode] cat6 ({tier_label}) start n={length(target_genes)}")

  class_df  <- internal_cache$classification
  mr_df     <- internal_cache$mr_results
  pqtl_pref <- internal_cache$pqtl_preferred

  rows <- lapply(target_genes, function(g) {
    class_row  <- if (!is.null(class_df) && "gene" %in% names(class_df))
      class_df[tolower(class_df$gene) == g, , drop = FALSE] else data.frame()
    gene_class <- failmode_get_chr(class_row, "class")

    pt_note         <- get_target_annotation(g, "post_transcriptional_note", diag_cfg)
    has_pt_mismatch <- !is.na(pt_note) && nzchar(pt_note)

    eqtl_not_supported <- !is.na(gene_class) && gene_class == "Not_supported"

    mr_primary <- if (!is.null(mr_df) && nrow(mr_df) > 0 && "gene" %in% names(mr_df)) {
      r <- mr_df[tolower(mr_df$gene) == g, , drop = FALSE]
      if ("outcome_name" %in% names(r)) {
        r2 <- r[r$outcome_name == "RA_EUR_PRIMARY", , drop = FALSE]
        if (nrow(r2) > 0) r2[1, , drop = FALSE] else if (nrow(r) > 0) r[1, , drop = FALSE] else data.frame()
      } else if (nrow(r) > 0) r[1, , drop = FALSE] else data.frame()
    } else data.frame()

    eqtl_beta <- failmode_get_num(mr_primary, "b")

    pqtl_primary <- if (!is.null(pqtl_pref) && nrow(pqtl_pref) > 0 && "gene" %in% names(pqtl_pref)) {
      r <- pqtl_pref[tolower(pqtl_pref$gene) == g, , drop = FALSE]
      if ("outcome_name" %in% names(r)) {
        r2 <- r[r$outcome_name == "RA_EUR_PRIMARY", , drop = FALSE]
        if (nrow(r2) > 0) r2[1, , drop = FALSE] else if (nrow(r) > 0) r[1, , drop = FALSE] else data.frame()
      } else if (nrow(r) > 0) r[1, , drop = FALSE] else data.frame()
    } else data.frame()

    pqtl_pval <- failmode_get_num(pqtl_primary, "pval")
    pqtl_beta <- failmode_get_num(pqtl_primary, "b")
    pqtl_sig  <- !is.na(pqtl_pval) && is.finite(pqtl_pval) && pqtl_pval < 0.05
    dir_concordant <- is.na(eqtl_beta) || is.na(pqtl_beta) || sign(eqtl_beta) == sign(pqtl_beta)
    has_discordant <- eqtl_not_supported && pqtl_sig && isTRUE(dir_concordant)

    # eQTL Catalogue sQTL OR QTLbase2 sQTL (no blood eQTL required)
    has_sqtl <- FALSE
    if (!is.null(qtl_cache) && !is.null(qtl_cache[[g]])) {
      has_sqtl <- isTRUE(qtl_cache[[g]]$eqtl_cat_sqtl_hits) ||
                  isTRUE(qtl_cache[[g]]$qtlbase2_has_sqtl)
    }
    gtex_blood_eqtl_g <- gtex_blood_eqtl_for_cat6(g, qtl_cache)
    sqtl_flag <- has_sqtl && !gtex_blood_eqtl_g

    # GoDMC cis-mQTL (mqtl_status == "tested_positive") OR QTLbase2 mQTL/hQTL/caQTL.
    has_alt_qtl <- FALSE
    if (!is.null(qtl_cache) && !is.null(qtl_cache[[g]])) {
      qc_g          <- qtl_cache[[g]]
      mqtl_status_g <- qc_g$mqtl_status %||% "data_unavailable"
      has_alt_qtl   <- identical(mqtl_status_g, "tested_positive") ||
                       isTRUE(qc_g$qtlbase2_has_mqtl)              ||
                       isTRUE(qc_g$qtlbase2_has_hqtl)              ||
                       isTRUE(qc_g$qtlbase2_has_caqtl)
    } else {
      mqtl_status_g <- "data_unavailable"
    }
    mqtl_flag <- has_alt_qtl && !gtex_blood_eqtl_g

    sqtl_discordance_flag <- gtex_blood_eqtl_g && has_sqtl

    # QTLbase2 xQTL types found (for output column)
    qtlbase2_xqtl_str <- if (!is.null(qtl_cache) && !is.null(qtl_cache[[g]])) {
      types <- qtl_cache[[g]]$qtlbase2_xqtl_types %||% character(0)
      if (length(types) > 0) paste(types, collapse = "; ") else NA_character_
    } else {
      NA_character_
    }

    # ── Priority: DISCORDANT > SQTL_DISCORDANCE > SPLICING > ALTERNATIVE > PT ──
    cat6_flag <- if (has_discordant)           "EQTL_PQTL_DISCORDANT"
    else if (sqtl_discordance_flag)            "CAT6_SQTL_DISCORDANCE"
    else if (sqtl_flag)                        "SPLICING_QTL_DETECTED"
    else if (mqtl_flag)                        "ALTERNATIVE_QTL_AVAILABLE"
    else if (has_pt_mismatch)                  "POST_TRANSCRIPTIONAL_MISMATCH"
    else                                       NA_character_

    cat6_confidence <- if (identical(cat6_flag, "EQTL_PQTL_DISCORDANT"))      "High"
    else if (identical(cat6_flag, "CAT6_SQTL_DISCORDANCE"))                   "Moderate"
    else if (identical(cat6_flag, "SPLICING_QTL_DETECTED"))                   "Moderate"
    else if (identical(cat6_flag, "ALTERNATIVE_QTL_AVAILABLE"))               "Moderate"
    else if (identical(cat6_flag, "POST_TRANSCRIPTIONAL_MISMATCH"))           "High"
    else                                                                       NA_character_

    alt_strategy <- lookup_alternative_strategy(cat6_flag, diag_cfg)

    # ── Implementation note: upgrade to "complete:" when qtl_cache available ──
    impl_note <- if (is.null(qtl_cache)) {
    } else if (has_sqtl) {
      "complete: sQTL confirmed"
    } else if (has_alt_qtl) {
      "complete: mQTL confirmed"
    } else if (mqtl_status_g %in% c("data_unavailable", "api_error")) {
      "complete: sQTL absent; mQTL data_unavailable (ARIES load failed)"
    } else {
      "complete: no sQTL/mQTL detected"
    }

    data.frame(
      gene                           = g,
      class                          = gene_class,
      cat6_flag                      = cat6_flag,
      cat6_confidence                = cat6_confidence,
      cat6_post_transcriptional_note = pt_note,
      cat6_eqtl_pqtl_discordant      = has_discordant,
      cat6_pqtl_primary_pval         = pqtl_pval,
      cat6_pqtl_primary_beta         = pqtl_beta,
      cat6_sqtl_detected             = sqtl_flag,
      cat6_sqtl_discordance          = sqtl_discordance_flag,
      cat6_mqtl_available            = mqtl_flag,
      cat6_qtlbase2_xqtl_types       = qtlbase2_xqtl_str,
      cat6_implementation_note       = impl_note,
      recommended_strategy           = alt_strategy,
      annotation_timestamp           = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
      stringsAsFactors               = FALSE
    )
  })

  result <- dplyr::bind_rows(rows)
  n_flagged <- sum(!is.na(result$cat6_flag), na.rm = TRUE)
  failmode_log(log_file,
    "[fail_mode] cat6 done n_genes={nrow(result)} n_flagged={n_flagged}")
  result
}

run_cat6_post_transcriptional_partial <- run_cat6_post_transcriptional

gtex_blood_eqtl_for_cat6 <- function(gene, qtl_cache) {
  if (is.null(qtl_cache) || is.null(qtl_cache[[gene]])) return(NA)
  isTRUE(qtl_cache[[gene]]$gtex_blood_eqtl)
}


# ─────────────────────────────────────────────────────────────────────────────
# COMPOSITE TABLE BUILDER
#
# cat_results: named list — cat7=df, cat3=df, cat1=df, cat5=df, cat6=df
# ─────────────────────────────────────────────────────────────────────────────

build_fail_mode_composite <- function(cat_results, diag_cfg) {
  cat7 <- cat_results[["cat7"]] %||% data.frame()
  if (is.null(cat7) || nrow(cat7) == 0) return(data.frame())

  out_cols <- diag_cfg$output$composite_table_columns %||%
    c("gene", "class",
      "cat1_flag", "cat1_confidence",
      "cat2_flag", "cat2_confidence", "cat2_locus_score",
      "cat3_flag", "cat3_confidence",
      "cat4_flag", "cat4_confidence",
      "cat5_flag", "cat5_confidence",
      "cat6_flag", "cat6_confidence",
      "cat7_flag", "cat7_confidence",
      "n_categories_flagged", "recommended_strategy", "annotation_timestamp")

  base_cols <- c("gene", "class", "cat7_flag", "cat7_confidence",
                 "recommended_strategy", "annotation_timestamp")
  composite <- cat7[, intersect(base_cols, names(cat7)), drop = FALSE]

  other_cats <- setdiff(names(cat_results), "cat7")
  for (cat_nm in other_cats) {
    cat_df <- cat_results[[cat_nm]]
    if (is.null(cat_df) || !is.data.frame(cat_df) || nrow(cat_df) == 0) next
    flag_col  <- paste0(cat_nm, "_flag")
    conf_col  <- paste0(cat_nm, "_confidence")
    extra_col <- if (cat_nm == "cat2") "cat2_locus_score" else character(0)
    join_cols <- c("gene", intersect(c(flag_col, conf_col, extra_col), names(cat_df)))
    if (length(join_cols) <= 1) next
    cat_sub <- cat_df[, join_cols, drop = FALSE]
    composite <- dplyr::left_join(composite, cat_sub, by = "gene")
  }

  for (col in setdiff(out_cols, names(composite))) {
    composite[[col]] <- NA
  }

  flag_cols <- grep("^cat[0-9]+_flag$", names(composite), value = TRUE)
  if (length(flag_cols) > 0) {
    composite$n_categories_flagged <- rowSums(
      !is.na(composite[, flag_cols, drop = FALSE]), na.rm = TRUE
    )
  }

  available_cols <- intersect(out_cols, names(composite))
  composite[, available_cols, drop = FALSE]
}


# ─────────────────────────────────────────────────────────────────────────────
# ORCHESTRATOR: run_failure_mode_diagnostic
#
#
# ─────────────────────────────────────────────────────────────────────────────

run_failure_mode_diagnostic <- function(results_dir,
                                        config_path,
                                        output_dir    = results_dir,
                                        log_file      = NULL,
                                        err_file      = NULL,
                                        run_api_fetch = FALSE,
                                        target_genes  = NULL) {
  failmode_log(log_file,
    "[fail_mode] START config={config_path} api_fetch={run_api_fetch}")

  diag_cfg <- tryCatch(
    load_fail_mode_config(config_path),
    error = function(e) {
      emsg <- conditionMessage(e)
      failmode_log(err_file %||% log_file, "[fail_mode] ERROR loading config: {emsg}")
      NULL
    }
  )
  if (is.null(diag_cfg)) {
    failmode_log(log_file, "[fail_mode] ABORTED config load failed")
    return(invisible(NULL))
  }

  internal_cache <- fetch_pipeline_internal(results_dir, log_file = log_file)
  if (!is.null(target_genes) && length(target_genes) > 0) {
    target_genes <- tolower(trimws(target_genes))
    failmode_log(log_file, "[fail_mode] target_genes override n={length(target_genes)} genes={paste(target_genes, collapse=';')}")
  } else {
    target_genes <- identify_diagnostic_targets(internal_cache, diag_cfg)
  }

  genomic_cache      <- NULL
  pharmacology_cache <- NULL
  qtl_cache          <- NULL

  api_log <- new_api_log()

  if (isTRUE(run_api_fetch) && length(target_genes) > 0) {
    genomic_cache <- tryCatch(
      fetch_genomic_context(target_genes, diag_cfg, log_file = log_file,
                            api_log = api_log),
      error = function(e) {
        emsg <- conditionMessage(e)
        failmode_log(err_file %||% log_file, "[fail_mode] ERROR fetch_genomic_context: {emsg}")
        NULL
      }
    )

    pharmacology_cache <- tryCatch(
      fetch_pharmacology(target_genes, diag_cfg, log_file = log_file,
                         genomic_cache = genomic_cache, api_log = api_log),
      error = function(e) {
        emsg <- conditionMessage(e)
        failmode_log(err_file %||% log_file, "[fail_mode] ERROR fetch_pharmacology: {emsg}")
        NULL
      }
    )

    # D: GTEx + eQTL Catalogue + GoDMC (cat4 + cat6 full)
    qtl_cache <- tryCatch(
      fetch_multiomics_qtl(target_genes, diag_cfg, log_file = log_file,
                           genomic_cache = genomic_cache, api_log = api_log,
                           lead_variants = internal_cache$lead_variants),
      error = function(e) {
        emsg <- conditionMessage(e)
        failmode_log(err_file %||% log_file, "[fail_mode] ERROR fetch_multiomics_qtl: {emsg}")
        NULL
      }
    )

    failmode_api_log_append(
      api_log,
      gene        = paste(target_genes, collapse = ";"),
      api_name    = "mqtl",
      endpoint    = "aries_mrinstruments_summary",
      cache_hit   = FALSE,
      success     = TRUE,
      http_status = NA_integer_,
      note        = "ARIES mQTL (MRInstruments) query completed"
    )
  }

  safe_run <- function(fn, label, ...) {
    tryCatch(fn(...), error = function(e) {
      emsg <- conditionMessage(e)
      failmode_log(err_file %||% log_file, paste0("[fail_mode] ERROR ", label, ": {emsg}"))
      data.frame()
    })
  }

  cat7_detail <- safe_run(run_cat7_iv_absent,          "cat7",
    internal_cache, diag_cfg, log_file = log_file)
  cat3_detail <- safe_run(run_cat3_pqtl_panel_absent,  "cat3",
    internal_cache, diag_cfg, log_file = log_file)
  cat1_detail <- safe_run(run_cat1_mhc_ld_complexity,  "cat1",
    internal_cache, diag_cfg, log_file = log_file)
  cat5_detail <- safe_run(run_cat5_mechanism_mismatch, "cat5",
    internal_cache, diag_cfg, log_file = log_file,
    pharmacology_cache = pharmacology_cache)
  cat6_detail <- safe_run(run_cat6_post_transcriptional, "cat6",
    internal_cache, diag_cfg, log_file = log_file,
    qtl_cache = qtl_cache)

  cat2_detail <- if (!is.null(genomic_cache)) {
    safe_run(run_cat2_structural_variant,      "cat2",
      internal_cache, genomic_cache, diag_cfg, log_file = log_file)
  } else data.frame()

  cat4_detail <- if (!is.null(qtl_cache)) {
    safe_run(run_cat4_tissue_cell_specificity, "cat4",
      internal_cache, qtl_cache, diag_cfg, log_file = log_file)
  } else data.frame()

  # ── Composite table ───────────────────────────────────────────────────────
  cat_results <- list(
    cat7 = cat7_detail, cat3 = cat3_detail, cat1 = cat1_detail,
    cat2 = cat2_detail, cat5 = cat5_detail, cat4 = cat4_detail,
    cat6 = cat6_detail
  )
  composite <- build_fail_mode_composite(cat_results, diag_cfg)

  out_cfg     <- diag_cfg$output$files %||% diag_cfg$output %||% list()
  f_composite <- as.character(out_cfg$composite_table %||% "fail_diagnostic_composite.csv")
  f_cat7    <- as.character(out_cfg$cat7_detail    %||% "cat7_iv_absent_detail.csv")
  f_cat3    <- as.character(out_cfg$cat3_detail    %||% "cat3_pqtl_panel_detail.csv")
  f_cat1    <- as.character(out_cfg$cat1_detail    %||% "cat1_mhc_ld_detail.csv")
  f_cat2    <- as.character(out_cfg$cat2_detail    %||% "cat2_paralog_ld_detail.csv")
  f_cat5    <- as.character(out_cfg$cat5_detail    %||% "cat5_mechanism_mismatch_detail.csv")
  f_cat4    <- as.character(out_cfg$cat4_detail    %||% "cat4_multiomics_detail.csv")
  f_cat6    <- as.character(out_cfg$cat6_detail    %||% "cat6_post_transcriptional_detail.csv")
  f_api_log <- as.character(out_cfg$api_call_log   %||% "fail_mode_annotation_log.csv")

  write_csv_safe(composite,   fs::path(output_dir, f_composite))
  write_csv_safe(cat7_detail, fs::path(output_dir, f_cat7))
  write_csv_safe(cat3_detail, fs::path(output_dir, f_cat3))
  write_csv_safe(cat1_detail, fs::path(output_dir, f_cat1))
  write_csv_safe(cat5_detail, fs::path(output_dir, f_cat5))
  write_csv_safe(cat6_detail, fs::path(output_dir, f_cat6))
  if (nrow(cat2_detail) > 0) write_csv_safe(cat2_detail, fs::path(output_dir, f_cat2))
  if (nrow(cat4_detail) > 0) write_csv_safe(cat4_detail, fs::path(output_dir, f_cat4))

  api_log_df <- failmode_api_log_to_df(api_log)
  if (nrow(api_log_df) > 0) write_csv_safe(api_log_df, fs::path(output_dir, f_api_log))

  failmode_log(log_file,
    "[fail_mode] DONE composite={f_composite} n_genes={nrow(composite)} api_fetch={run_api_fetch} api_calls={nrow(api_log_df)}")

  invisible(list(
    composite = composite,
    cat7 = cat7_detail, cat3 = cat3_detail, cat1 = cat1_detail,
    cat2 = cat2_detail, cat5 = cat5_detail, cat4 = cat4_detail,
    cat6 = cat6_detail,
    api_log = api_log_df
  ))
}


# ══════════════════════════════════════════════════════════════════════════════
#                         TIER B — EXTERNAL API LAYER
# ══════════════════════════════════════════════════════════════════════════════

# ─────────────────────────────────────────────────────────────────────────────
# API UTILITY HELPERS
#
# ─────────────────────────────────────────────────────────────────────────────

failmode_api_cfg <- function(diag_cfg, delay_override_ms = NULL) {
  gl  <- diag_cfg$api_global %||% list()
  rl  <- gl$rate_limiting   %||% list()
  ret <- gl$retry           %||% list()
  htp <- gl$http            %||% list()
  list(
    delay_ms       = delay_override_ms %||% rl$default_delay_ms %||% 200,
    max_retries    = ret$max_retries   %||% 3L,
    timeout_sec    = htp$timeout_sec   %||% 30L,
    retry_delay_s  = ret$retry_delay_sec  %||% 5L,
    retryable      = ret$retryable_http_codes %||% c(429L, 500L, 502L, 503L, 504L),
    user_agent     = htp$user_agent    %||% "btsDMARDs_MR_pipeline/v2.0"
  )
}

failmode_cache_path <- function(api_name, gene, diag_cfg) {
  cache_dir <- diag_cfg$api_global$cache$directory %||% "cache/failmode_api"
  fs::path(cache_dir, paste0(api_name, "_", gene, ".rds"))
}

failmode_cache_read <- function(api_name, gene, diag_cfg) {
  if (!isTRUE(diag_cfg$api_global$cache$enabled %||% TRUE)) return(NULL)
  path <- failmode_cache_path(api_name, gene, diag_cfg)
  if (!fs::file_exists(path)) return(NULL)
  expiry_days <- diag_cfg$api_global$cache$expiry_days %||% 30
  age_days <- as.numeric(difftime(Sys.time(), fs::file_info(path)$modification_time, units = "days"))
  if (!is.na(age_days) && age_days > expiry_days) return(NULL)
  tryCatch(readRDS(path), error = function(e) NULL)
}

failmode_cache_write <- function(obj, api_name, gene, diag_cfg) {
  if (!isTRUE(diag_cfg$api_global$cache$enabled %||% TRUE)) return(invisible(NULL))
  path <- failmode_cache_path(api_name, gene, diag_cfg)
  tryCatch({
    fs::dir_create(fs::path_dir(path), recurse = TRUE)
    saveRDS(obj, path)
  }, error = function(e) invisible(NULL))
}

# HTTP GET with rate-limit + retry
failmode_api_get <- function(url, cfg) {
  Sys.sleep(cfg$delay_ms / 1000)
  for (attempt in seq_len(cfg$max_retries)) {
    resp <- tryCatch(
      httr::GET(url,
                httr::timeout(cfg$timeout_sec),
                httr::user_agent(cfg$user_agent),
                httr::accept_json()),
      error = function(e) NULL
    )
    if (is.null(resp)) break
    sc <- httr::status_code(resp)
    if (sc == 200L) {
      return(tryCatch(
        httr::content(resp, as = "parsed", type = "application/json", encoding = "UTF-8"),
        error = function(e) NULL
      ))
    }
    if (!sc %in% cfg$retryable) break
    if (attempt < cfg$max_retries) Sys.sleep(cfg$retry_delay_s)
  }
  NULL
}

# ── API Call Log (reference-semantics accumulator) ────────────────────────

new_api_log <- function() {
  e <- new.env(parent = emptyenv())
  e$entries <- list()
  e
}

failmode_api_log_append <- function(api_log, gene, api_name, endpoint,
                                     cache_hit, http_status = NA_integer_,
                                     duration_ms = NA_real_, success,
                                     note = NA_character_) {
  if (is.null(api_log)) return(invisible(NULL))
  api_log$entries <- c(api_log$entries, list(data.frame(
    timestamp    = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    gene         = as.character(gene),
    api_name     = as.character(api_name),
    endpoint     = as.character(endpoint),
    cache_hit    = as.logical(cache_hit),
    http_status  = suppressWarnings(as.integer(http_status)),
    duration_ms  = as.numeric(duration_ms),
    success      = as.logical(success),
    note         = as.character(note),
    stringsAsFactors = FALSE
  )))
  invisible(NULL)
}

failmode_api_log_to_df <- function(api_log) {
  empty <- data.frame(
    timestamp = character(0), gene = character(0), api_name = character(0),
    endpoint = character(0), cache_hit = logical(0), http_status = integer(0),
    duration_ms = numeric(0), success = logical(0), note = character(0),
    stringsAsFactors = FALSE
  )
  if (is.null(api_log) || length(api_log$entries) == 0) return(empty)
  dplyr::bind_rows(api_log$entries)
}


failmode_timed_get <- function(url, cfg, api_name, gene, api_log) {
  t0 <- proc.time()[["elapsed"]]
  result <- failmode_api_get(url, cfg)
  dur <- (proc.time()[["elapsed"]] - t0) * 1000
  failmode_api_log_append(api_log, gene, api_name, url,
    cache_hit = FALSE, duration_ms = dur, success = !is.null(result))
  result
}

failmode_timed_graphql <- function(url, query_str, variables = list(), cfg,
                                    api_name, gene, api_log) {
  t0 <- proc.time()[["elapsed"]]
  result <- failmode_api_graphql(url, query_str, variables, cfg)
  dur <- (proc.time()[["elapsed"]] - t0) * 1000
  failmode_api_log_append(api_log, gene, api_name,
    paste0(url, " [GraphQL]"),
    cache_hit = FALSE, duration_ms = dur, success = !is.null(result))
  result
}

# GraphQL POST (Open Targets, gnomAD)
failmode_api_graphql <- function(url, query_str, variables = list(), cfg) {
  body <- jsonlite::toJSON(list(query = query_str, variables = variables), auto_unbox = TRUE)
  Sys.sleep(cfg$delay_ms / 1000)
  for (attempt in seq_len(cfg$max_retries)) {
    resp <- tryCatch(
      httr::POST(url,
                 httr::timeout(cfg$timeout_sec),
                 httr::user_agent(cfg$user_agent),
                 httr::content_type_json(),
                 httr::accept_json(),
                 body = body),
      error = function(e) NULL
    )
    if (is.null(resp)) break
    sc <- httr::status_code(resp)
    if (sc == 200L) {
      return(tryCatch(
        httr::content(resp, as = "parsed", type = "application/json", encoding = "UTF-8"),
        error = function(e) NULL
      ))
    }
    if (!sc %in% cfg$retryable) break
    if (attempt < cfg$max_retries) Sys.sleep(cfg$retry_delay_s)
  }
  NULL
}


# ─────────────────────────────────────────────────────────────────────────────
# FETCH GROUP B: GENOMIC CONTEXT
#
# Ensembl GRCh37 REST → paralog → gnomAD SV GraphQL → UCSC mappability
# ─────────────────────────────────────────────────────────────────────────────

fetch_genomic_context <- function(target_genes, diag_cfg, log_file = NULL,
                                  api_log = NULL) {
  failmode_log(log_file,
    "[fail_mode] fetch_genomic_context start n={length(target_genes)}")

  ens_cfg  <- failmode_api_cfg(diag_cfg,
    delay_override_ms = diag_cfg$api_global$rate_limiting$ensembl_delay_ms %||% 350)
  rest_cfg <- failmode_api_cfg(diag_cfg)
  mhc_core <- diag_cfg$fetch_genomic_context$mhc_boundaries$grch37$core %||%
    list(chr = "6", start = 28477797L, end = 33448354L)

  # ── Panel coordinate lookup table (Option A: cross-reference for n_paralogs_nearby) ──
  # Ensembl /homology/symbol/ returns target$id (Ensembl gene ID) but NOT genomic
  # coordinates. Pre-load all panel genes' ensembl_gene cache entries, keyed by
  # Ensembl ID, so we can resolve coordinates for in-panel paralogs without extra
  # API calls.
  panel_coords <- local({
    pc <- list()
    for (pg in target_genes) {
      ce <- failmode_cache_read("ensembl_gene", pg, diag_cfg)
      if (!is.null(ce) && nzchar(ce$id %||% "")) {
        pc[[ce$id]] <- list(
          gene  = pg,
          chr   = gsub("^chr", "", as.character(ce$seq_region_name %||% "")),
          start = suppressWarnings(as.integer(ce$start %||% NA_integer_))
        )
      }
    }
    pc
  })

  result <- lapply(target_genes, function(g) {
    gene_upper <- toupper(g)

    # ── Step 1: Ensembl gene lookup (GRCh37) ──────────────────────────────
    ens_entry <- failmode_cache_read("ensembl_gene", g, diag_cfg)
    if (!is.null(ens_entry)) {
      failmode_api_log_append(api_log, g, "ensembl_gene",
        paste0("https://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens/", toupper(g)),
        cache_hit = TRUE, success = TRUE)
    } else {
      url <- paste0("https://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens/",
                    gene_upper, "?content-type=application/json&expand=0")
      ens_entry <- failmode_timed_get(url, ens_cfg, "ensembl_gene", g, api_log)
      if (!is.null(ens_entry)) failmode_cache_write(ens_entry, "ensembl_gene", g, diag_cfg)
    }

    ensembl_id <- as.character(ens_entry$id                %||% NA_character_)
    chr        <- as.character(ens_entry$seq_region_name   %||% NA_character_)
    start      <- suppressWarnings(as.integer(ens_entry$start %||% NA_integer_))
    end        <- suppressWarnings(as.integer(ens_entry$end   %||% NA_integer_))

    # ── Step 2: Paralog query ─────────────────────────────────────────────
    n_paralogs_nearby    <- NA_integer_
    n_paralogs_all       <- NA_integer_
    paralogs_nearby      <- character(0)
    max_paralog_identity <- NA_real_

    # Paralog query uses gene symbol (grch37 homology/id returns 404; symbol path reliable)
    if (TRUE) {
      par_entry <- failmode_cache_read("ensembl_paralog", g, diag_cfg)
      par_url   <- paste0("https://rest.ensembl.org/homology/symbol/homo_sapiens/", gene_upper,
                          "?type=paralogues&content-type=application/json")
      if (!is.null(par_entry)) {
        failmode_api_log_append(api_log, g, "ensembl_paralog", par_url,
          cache_hit = TRUE, success = TRUE)
      } else {
        par_entry <- failmode_timed_get(par_url, ens_cfg, "ensembl_paralog", g, api_log)
        if (!is.null(par_entry)) failmode_cache_write(par_entry, "ensembl_paralog", g, diag_cfg)
      }
      if (!is.null(par_entry$data) && length(par_entry$data) > 0) {
        homologies <- par_entry$data[[1]]$homologies %||% list()
        n_paralogs_all <- length(homologies)
        # Extract max sequence identity across all paralogs
        perc_ids <- vapply(homologies, function(h) {
          suppressWarnings(as.numeric(h$target$perc_id %||% NA_real_))
        }, numeric(1))
        if (any(!is.na(perc_ids))) max_paralog_identity <- max(perc_ids, na.rm = TRUE)
        # Nearby paralogs (within 1 Mb on same chromosome) — Option A
        # Ensembl homology API does not include coordinates in target fields.
        # Strategy: look up each paralog's Ensembl ID in panel_coords (pre-built
        # above from ensembl_gene cache); fall back to API fields (usually absent)
        # for non-panel paralogs.
        g_chr_norm <- gsub("^chr", "", as.character(chr %||% ""))
        nearby_flags <- vapply(homologies, function(h) {
          tgt_id <- as.character(h$target$id %||% "")
          pc_hit <- if (nzchar(tgt_id)) panel_coords[[tgt_id]] else NULL
          if (!is.null(pc_hit)) {
            t_chr   <- pc_hit$chr
            t_start <- pc_hit$start
          } else {
            # Fallback: API-returned coords (typically absent for non-panel paralogs)
            tgt     <- h$target %||% list()
            t_chr   <- gsub("^chr", "", as.character(tgt$seq_region_name %||% ""))
            t_start <- suppressWarnings(as.integer(tgt$start %||% NA_integer_))
          }
          !is.na(start) && !is.na(t_start) && nzchar(g_chr_norm) && nzchar(t_chr) &&
            t_chr == g_chr_norm && abs(t_start - start) <= 1000000L
        }, logical(1))
        n_paralogs_nearby <- sum(nearby_flags, na.rm = TRUE)
        # Report gene symbol for panel paralogs, Ensembl ID otherwise
        paralogs_nearby <- vapply(homologies[nearby_flags], function(h) {
          tgt_id <- as.character(h$target$id %||% NA_character_)
          pc_hit <- if (!is.na(tgt_id) && nzchar(tgt_id)) panel_coords[[tgt_id]] else NULL
          if (!is.null(pc_hit)) pc_hit$gene else tgt_id
        }, character(1))
        paralogs_nearby <- paralogs_nearby[!is.na(paralogs_nearby)]
      }
    }

    # ── Step 3: gnomAD SV (GraphQL) ───────────────────────────────────────
    sv_count <- NA_integer_

    if (!is.na(chr) && !is.na(start) && !is.na(end)) {
      sv_entry <- failmode_cache_read("gnomad_sv", g, diag_cfg)
      if (!is.null(sv_entry)) {
        failmode_api_log_append(api_log, g, "gnomad_sv",
          "https://gnomad.broadinstitute.org/api [GraphQL]",
          cache_hit = TRUE, success = TRUE)
      } else {
        # gnomAD v4 API: region requires reference_genome; structural_variants requires dataset arg
        # GRCh37 coords used as approximation (SVs are large; coordinate shift is negligible)
        gql <- sprintf(
          '{region(chrom:"%s",start:%d,stop:%d,reference_genome:GRCh38){structural_variants(dataset:gnomad_sv_r4){variant_id type}}}',
          chr, max(1L, start - 500000L), end + 500000L
        )
        sv_entry <- failmode_timed_graphql(
          "https://gnomad.broadinstitute.org/api", gql, list(), rest_cfg,
          "gnomad_sv", g, api_log)
        if (!is.null(sv_entry)) failmode_cache_write(sv_entry, "gnomad_sv", g, diag_cfg)
      }
      svs <- sv_entry$data$region$structural_variants %||% list()
      sv_count <- length(svs)
    }

    # ── Step 4: UCSC mappability ──────────────────────────────────────────
    mappability_mean <- NA_real_

    if (!is.na(chr) && !is.na(start) && !is.na(end)) {
      map_entry <- failmode_cache_read("ucsc_map", g, diag_cfg)
      if (!is.null(map_entry)) {
        failmode_api_log_append(api_log, g, "ucsc_map",
          sprintf("https://api.genome.ucsc.edu/getData/track?genome=hg19&track=wgEncodeCrgMapabilityAlign36mer&chrom=chr%s", chr),
          cache_hit = TRUE, success = TRUE)
      } else {
        url3 <- sprintf(
          "https://api.genome.ucsc.edu/getData/track?genome=hg19&track=wgEncodeCrgMapabilityAlign36mer&chrom=chr%s&start=%d&end=%d",
          chr, start, end)
        map_entry <- failmode_timed_get(url3, rest_cfg, "ucsc_map", g, api_log)
        if (!is.null(map_entry)) failmode_cache_write(map_entry, "ucsc_map", g, diag_cfg)
      }
      scores <- tryCatch({
        raw <- map_entry$wgEncodeCrgMapabilityAlign36mer %||% list()
        vals <- vapply(raw, function(x) suppressWarnings(as.numeric(x$score %||% NA_real_)), numeric(1))
        vals[!is.na(vals)]
      }, error = function(e) numeric(0))
      if (length(scores) > 0) mappability_mean <- mean(scores)
    }

    list(gene = g, ensembl_id = ensembl_id, chr = chr, start = start, end = end,
         n_paralogs_nearby = n_paralogs_nearby, n_paralogs_all = n_paralogs_all,
         paralogs_nearby = paralogs_nearby, max_paralog_identity = max_paralog_identity,
         sv_count = sv_count, mappability_mean = mappability_mean)
  })

  names(result) <- target_genes
  n_ok <- sum(vapply(result, function(x) !is.na(x$ensembl_id %||% NA), logical(1)))
  failmode_log(log_file,
    "[fail_mode] fetch_genomic_context done ensembl_resolved={n_ok}/{length(target_genes)}")
  result
}


# ─────────────────────────────────────────────────────────────────────────────
# FETCH GROUP C: PHARMACOLOGY
#
# UniProt REST → Open Targets GraphQL
# ─────────────────────────────────────────────────────────────────────────────

fetch_pharmacology <- function(target_genes, diag_cfg, log_file = NULL,
                               genomic_cache = NULL, api_log = NULL) {
  failmode_log(log_file,
    "[fail_mode] fetch_pharmacology start n={length(target_genes)}")

  rest_cfg <- failmode_api_cfg(diag_cfg)
  ot_cfg   <- failmode_api_cfg(diag_cfg,
    delay_override_ms = diag_cfg$api_global$rate_limiting$open_targets_delay_ms %||% 100)

  result <- lapply(target_genes, function(g) {
    uniprot_id <- get_target_annotation(g, "uniprot_id", diag_cfg)
    ensembl_id <- if (!is.null(genomic_cache[[g]])) genomic_cache[[g]]$ensembl_id else NA_character_

    # ── UniProt: subcellular location ─────────────────────────────────────
    uniprot_subcell <- NA_character_
    uniprot_go_mf   <- NA_character_

    if (!is.na(uniprot_id)) {
      up_entry <- failmode_cache_read("uniprot", g, diag_cfg)
      if (!is.null(up_entry)) {
        failmode_api_log_append(api_log, g, "uniprot",
          paste0("https://rest.uniprot.org/uniprotkb/", uniprot_id),
          cache_hit = TRUE, success = TRUE)
      } else {
        url <- paste0("https://rest.uniprot.org/uniprotkb/", uniprot_id,
                      "?fields=cc_subcellular_location,keyword&format=json")
        up_entry <- failmode_timed_get(url, rest_cfg, "uniprot", g, api_log)
        if (!is.null(up_entry)) failmode_cache_write(up_entry, "uniprot", g, diag_cfg)
      }
      if (!is.null(up_entry)) {
        # Subcellular location
        comments <- up_entry$comments %||% list()
        subcell_c <- Filter(function(c) identical(c$commentType, "SUBCELLULAR LOCATION"), comments)
        if (length(subcell_c) > 0) {
          locs <- subcell_c[[1]]$subcellularLocations %||% list()
          loc_names <- vapply(locs, function(l) as.character(l$location$value %||% ""), character(1))
          uniprot_subcell <- paste(loc_names[nzchar(loc_names)], collapse = "; ")
          if (!nzchar(uniprot_subcell)) uniprot_subcell <- NA_character_
        }
        # GO molecular function (from keywords)
        kw <- up_entry$keywords %||% list()
        mf_kw <- Filter(function(k) identical(k$category, "Molecular function"), kw)
        if (length(mf_kw) > 0) {
          mf_names <- vapply(mf_kw, function(k) as.character(k$name %||% ""), character(1))
          uniprot_go_mf <- paste(mf_names[nzchar(mf_names)], collapse = "; ")
          if (!nzchar(uniprot_go_mf)) uniprot_go_mf <- NA_character_
        }
      }
    }

    # ── Open Targets: drug mechanism + tractability ───────────────────────
    ot_mechanisms    <- NA_character_
    ot_tractability  <- NA_character_
    ot_highest_phase <- NA_integer_

    if (!is.na(ensembl_id)) {
      ot_entry <- failmode_cache_read("open_targets", g, diag_cfg)
      if (!is.null(ot_entry)) {
        failmode_api_log_append(api_log, g, "open_targets",
          "https://api.platform.opentargets.org/api/v4/graphql [GraphQL]",
          cache_hit = TRUE, success = TRUE)
      } else {
        # OT Platform v24+ API: knownDrugs → drugAndClinicalCandidates,
        #   maximumClinicalTrialPhase → maxClinicalStage (string enum),
        #   tractability{id} → tractability{label}
        ot_query <- 'query($id:String!){target(ensemblId:$id){
          drugAndClinicalCandidates{rows{maxClinicalStage
            drug{mechanismsOfAction{rows{actionType mechanismOfAction}}}}}
          tractability{modality label value}
          subcellularLocations{location}}}'
        ot_entry <- failmode_timed_graphql(
          "https://api.platform.opentargets.org/api/v4/graphql",
          ot_query, list(id = ensembl_id), ot_cfg,
          "open_targets", g, api_log)
        if (!is.null(ot_entry)) failmode_cache_write(ot_entry, "open_targets", g, diag_cfg)
      }
      if (!is.null(ot_entry$data$target)) {
        drug_rows <- ot_entry$data$target$drugAndClinicalCandidates$rows %||% list()
        if (length(drug_rows) > 0) {
          # maxClinicalStage: string enum (APPROVAL, PHASE_4..PHASE_0, PHASE_2_3, PHASE_1_2, etc.)
          # Use ceiling of range for combined phases (e.g. PHASE_2_3 → 3)
          stage_to_int <- c(APPROVAL = 4L, PHASE_4 = 4L,
                            PHASE_3 = 3L, PHASE_2_3 = 3L,
                            PHASE_2 = 2L, PHASE_1_2 = 2L,
                            PHASE_1 = 1L, PHASE_0 = 0L)
          phases <- vapply(drug_rows, function(d) {
            s <- as.character(d$maxClinicalStage %||% "")
            stage_to_int[s][[1L]] %||% NA_integer_
          }, integer(1))
          if (any(!is.na(phases))) ot_highest_phase <- max(phases, na.rm = TRUE)
          action_types <- unique(unlist(lapply(drug_rows, function(d)
            vapply(d$drug$mechanismsOfAction$rows %||% list(),
                   function(m) as.character(m$actionType %||% ""), character(1)))))
          action_types <- action_types[nzchar(action_types)]
          ot_mechanisms <- if (length(action_types) > 0) paste(action_types, collapse = "; ") else NA_character_
        }
        tract_rows <- ot_entry$data$target$tractability %||% list()
        active_mod <- vapply(tract_rows, function(t)
          if (isTRUE(t$value)) as.character(t$modality %||% "") else "", character(1))
        active_mod <- unique(active_mod[nzchar(active_mod)])
        ot_tractability <- if (length(active_mod) > 0) paste(active_mod, collapse = "; ") else NA_character_
      }
    }

    list(gene = g, uniprot_id = uniprot_id,
         uniprot_subcellular_location = uniprot_subcell,
         uniprot_go_molecular_function = uniprot_go_mf,
         ot_mechanisms = ot_mechanisms,
         ot_tractability = ot_tractability,
         ot_highest_phase = ot_highest_phase)
  })

  names(result) <- target_genes
  n_ok <- sum(vapply(result, function(x) !is.na(x$ot_mechanisms %||% NA), logical(1)))
  failmode_log(log_file,
    "[fail_mode] fetch_pharmacology done ot_resolved={n_ok}/{length(target_genes)}")
  result
}


# ─────────────────────────────────────────────────────────────────────────────
# CHROMATIN STATE HELPERS (Roadmap/ENCODE ChromHMM BED)
#
# Pre-downloaded 15-state ChromHMM BED.gz files → TSS ± 2kb state lookup
# ─────────────────────────────────────────────────────────────────────────────

parse_chromhmm_label <- function(label) {
  if (is.null(label) || length(label) == 0 || is.na(label)) return(NA_character_)
  as.character(sub("^[0-9]+_", "", label[1]))
}

# Returns named list: blood/tcell/bcell/monocyte → data.frame(chr, start, end, state)
load_chromatin_beds <- function(diag_cfg, log_file = NULL) {
  chr_cfg  <- diag_cfg$fetch_multiomics_qtl$chromatin_query %||% list()
  bed_dir  <- chr_cfg$bed_directory %||% "data/roadmap_encode"
  bed_files <- chr_cfg$bed_files %||% list(
    blood    = "E062_15_coreMarks_mnemonics.bed.gz",
    tcell    = "E034_15_coreMarks_mnemonics.bed.gz",
    bcell    = "E045_15_coreMarks_mnemonics.bed.gz",
    monocyte = "E029_15_coreMarks_mnemonics.bed.gz"
  )

  result <- lapply(names(bed_files), function(tissue) {
    path <- file.path(bed_dir, bed_files[[tissue]])
    if (!file.exists(path)) {
      failmode_log(log_file,
        "[fail_mode] chromatin BED not found ({tissue}): {path}")
      return(NULL)
    }
    df <- tryCatch({
      con <- gzfile(path, open = "rt")
      on.exit(close(con), add = TRUE)
      read.table(con, sep = "\t", header = FALSE,
                 col.names = c("chr", "start", "end", "state"),
                 colClasses = c("character", "integer", "integer", "character"),
                 comment.char = "")
    }, error = function(e) {
      failmode_log(log_file,
        "[fail_mode] chromatin BED read error ({tissue}): {conditionMessage(e)}")
      NULL
    })
    if (!is.null(df)) {
      failmode_log(log_file,
        "[fail_mode] chromatin loaded {tissue} n_rows={nrow(df)}")
    }
    df
  })
  names(result) <- names(bed_files)

  if (all(vapply(result, is.null, logical(1)))) {
    failmode_log(log_file, "[fail_mode] chromatin: no BED files loaded")
    return(NULL)
  }
  result
}

# Returns raw label (e.g. "13_ReprPC") or NA
get_chromatin_state_at_tss <- function(gene_chr, gene_start, bed_df,
                                        window_bp = 2000L) {
  if (is.null(bed_df) || is.na(gene_chr) || is.na(gene_start)) return(NA_character_)
  tss         <- as.integer(gene_start)
  qstart      <- max(0L, tss - window_bp)
  qend        <- tss + window_bp
  bed_chr_fmt <- paste0("chr", gene_chr)     # Ensembl "6" → BED "chr6"

  rows <- bed_df[
    bed_df$chr == bed_chr_fmt &
    bed_df$end   > qstart    &
    bed_df$start < qend,
    , drop = FALSE
  ]
  if (nrow(rows) == 0L) return(NA_character_)

  rows$ovlp <- pmin(rows$end, qend) - pmax(rows$start, qstart)
  rows$state[which.max(rows$ovlp)]
}


# ─────────────────────────────────────────────────────────────────────────────
# QTLbase2 LOCAL CSV LOADER
#
# ─────────────────────────────────────────────────────────────────────────────

load_qtlbase2_local <- function(diag_cfg, log_file = NULL) {
  q2_cfg   <- diag_cfg$fetch_multiomics_qtl$qtlbase2 %||% list()
  enabled  <- isTRUE(q2_cfg$enabled %||% TRUE)
  csv_path <- as.character(q2_cfg$local_csv_path %||% "data/qtlbase2_lookup.csv")
  min_hits <- as.integer(q2_cfg$significance_filter$min_hits %||% 1L)

  if (!enabled) {
    failmode_log(log_file, "[fail_mode] qtlbase2 disabled, skipping local CSV")
    return(NULL)
  }

  if (!file.exists(csv_path)) {
    failmode_log(log_file,
      "[fail_mode] qtlbase2 local CSV not found: {csv_path} — proceeding without")
    return(NULL)
  }

  df <- tryCatch(
    utils::read.csv(csv_path, stringsAsFactors = FALSE, comment.char = "#",
                    na.strings = c("", "NA")),
    error = function(e) {
      failmode_log(log_file,
        "[fail_mode] qtlbase2 CSV read error: {conditionMessage(e)}")
      NULL
    }
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)

  # Validate required columns
  req_cols <- c("gene", "qtl_type", "n_hits")
  miss     <- setdiff(req_cols, names(df))
  if (length(miss) > 0) {
    failmode_log(log_file,
      "[fail_mode] qtlbase2 CSV missing columns: {paste(miss, collapse=', ')}")
    return(NULL)
  }

  df$gene     <- tolower(trimws(as.character(df$gene)))
  df$qtl_type <- trimws(as.character(df$qtl_type))
  df$n_hits   <- suppressWarnings(as.integer(df$n_hits))

  # Retain only rows with at least min_hits significant associations
  df <- df[!is.na(df$n_hits) & df$n_hits >= min_hits, , drop = FALSE]

  supported <- c("mQTL", "hQTL", "caQTL", "sQTL")
  df <- df[df$qtl_type %in% supported, , drop = FALSE]

  if (nrow(df) == 0) {
    failmode_log(log_file, "[fail_mode] qtlbase2 CSV loaded but no significant hits")
    return(NULL)
  }

  # Build named list per gene
  genes_in_csv <- unique(df$gene)
  result <- lapply(genes_in_csv, function(g) {
    sub <- df[df$gene == g, , drop = FALSE]
    types_found <- unique(sub$qtl_type)
    list(
      has_mqtl  = "mQTL"  %in% types_found,
      has_hqtl  = "hQTL"  %in% types_found,
      has_caqtl = "caQTL" %in% types_found,
      has_sqtl  = "sQTL"  %in% types_found,
      xqtl_types = types_found,
      records   = sub
    )
  })
  names(result) <- genes_in_csv

  n_genes  <- length(result)
  n_hits_t <- nrow(df)
  failmode_log(log_file,
    "[fail_mode] qtlbase2 loaded csv={csv_path} genes={n_genes} total_hits={n_hits_t}")
  result
}


# ─────────────────────────────────────────────────────────────────────────────
# ARIES mQTL HELPER  (MRInstruments v0.3.2)
#
#   - 5 developmental timepoints: Birth / 6mo / 17mo / 7yr / adolescence
#
#
# ─────────────────────────────────────────────────────────────────────────────

fetch_aries_mqtl <- function(gene_symbol, aries_data, pval_threshold = 1e-5) {
  # aries_data: pre-loaded data.frame (MRInstruments::aries_mqtl)
  # gene_symbol: lowercase (pipeline internal) — converted to UPPER for lookup
  gene_upper <- toupper(gene_symbol)

  hits <- aries_data[aries_data$gene == gene_upper &
                       !is.na(aries_data$pval) &
                       aries_data$pval < pval_threshold, ]

  n_hits <- nrow(hits)
  list(
    mqtl_status = if (n_hits > 0L) "tested_positive" else "tested_negative",
    n_hits      = n_hits,
    top_snp     = if (n_hits > 0L) hits$SNP[which.min(hits$pval)]  else NA_character_,
    top_cpg     = if (n_hits > 0L) hits$cpg[which.min(hits$pval)]  else NA_character_,
    top_p       = if (n_hits > 0L) min(hits$pval, na.rm = TRUE)    else NA_real_,
    source      = "aries_mrinstruments",
    note        = paste0(
        "ARIES blood mQTL (MRInstruments v",
        packageVersion("MRInstruments"),
        "); 5 timepoints (Birth-adolescence)"
    )
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# FETCH GROUP D: MULTI-OMICS QTL
#
# GTEx v2 → eQTL Catalogue v2 (eQTL + sQTL) → ieugwasr mQTL → QTLbase2 local CSV
# ─────────────────────────────────────────────────────────────────────────────

fetch_multiomics_qtl <- function(target_genes, diag_cfg, log_file = NULL,
                                 genomic_cache = NULL, api_log = NULL,
                                 lead_variants = NULL) {
  failmode_log(log_file,
    "[fail_mode] fetch_multiomics_qtl start n={length(target_genes)}")

  qtl_cfg   <- diag_cfg$fetch_multiomics_qtl %||% list()
  rest_cfg  <- failmode_api_cfg(diag_cfg)
  eqtl_cfg  <- failmode_api_cfg(diag_cfg,
    delay_override_ms = diag_cfg$api_global$rate_limiting$eqtl_catalogue_delay_ms %||% 500)

  gtex_tissues    <- qtl_cfg$gtex_query$iterate_over_datasets %||%
    qtl_cfg$gtex_query$iterate_over %||%
    c("Whole_Blood", "Spleen", "Cells_EBV-transformed_lymphocytes",
      "Cells_Cultured_fibroblasts", "Adipose_Subcutaneous")
  # eQTL Catalogue v2: QTD-format dataset IDs
  # ge (gene expression eQTL) — immune cell types
  eqtl_ge_datasets <- qtl_cfg$eqtl_catalogue_query$iterate_over_ge_datasets %||%
    c("QTD000021",   # BLUEPRINT monocyte ge
      "QTD000022",   # BLUEPRINT T-cell ge
      "QTD000023")   # BLUEPRINT neutrophil ge

  # sQTL (leafcutter / txrev) — immune cell types; confirmed in eQTL Catalogue v2
  eqtl_sqtl_datasets <- qtl_cfg$eqtl_catalogue_query$iterate_over_sqtl_datasets %||%
    c("QTD000359",   # GTEx v8 whole blood txrev (confirmed HTTP 200)
      "QTD000360",   # GTEx v8 whole blood leafcutter
      "QTD000025",   # BLUEPRINT monocyte leafcutter
      "QTD000035",   # BLUEPRINT CD4T leafcutter
      "QTD000478",   # Schmiedel 2018 B cell leafcutter
      "QTD000553")   # eQTL Catalogue additional immune sQTL

  # ── QTLbase2: local CSV (no API, loaded once for all genes) ───────────────
  qtlbase2_data <- load_qtlbase2_local(diag_cfg, log_file = log_file)

  aries_data <- tryCatch({
    if (!requireNamespace("MRInstruments", quietly = TRUE))
      stop("MRInstruments not installed")
    utils::data("aries_mqtl", package = "MRInstruments", envir = environment())
    get("aries_mqtl", envir = environment())
  }, error = function(e) {
    failmode_log(log_file,
      "[fail_mode] ARIES mQTL load failed: {conditionMessage(e)} — mqtl_status=data_unavailable")
    NULL
  })
  aries_cfg      <- qtl_cfg$aries_query %||% list()
  aries_pval_thr <- as.numeric(aries_cfg$pval_threshold %||% 1e-5)
  if (!is.null(aries_data)) {
    failmode_log(log_file,
      "[fail_mode] ARIES mQTL loaded: nrow={nrow(aries_data)} pval_threshold={aries_pval_thr}")
  }

  chr_cfg     <- qtl_cfg$chromatin_query %||% list()
  chr_enabled <- isTRUE(chr_cfg$enabled %||% TRUE)
  active_mks   <- chr_cfg$chromatin_states$active_marks    %||%
    c("TssA", "TssAFlnk", "EnhA1", "EnhA2")
  repressed_mks <- chr_cfg$chromatin_states$repressed_marks %||%
    c("ReprPC", "ReprPCWk")
  tss_window  <- as.integer(chr_cfg$tss_window_bp %||% 2000L)

  chromatin_beds <- if (chr_enabled) {
    needs_bed <- Filter(
      function(g) is.null(failmode_cache_read("chromatin", g, diag_cfg)),
      target_genes
    )
    if (length(needs_bed) > 0) {
      load_chromatin_beds(diag_cfg, log_file = log_file)
    } else {
      failmode_log(log_file, "[fail_mode] chromatin: all genes in cache, skipping BED load")
      NULL
    }
  } else NULL

  result <- lapply(target_genes, function(g) {
    gc         <- genomic_cache[[g]] %||% list()
    ensembl_id <- gc$ensembl_id %||% NA_character_
    gene_chr   <- gc$chr        %||% NA_character_
    gene_start <- gc$start      %||% NA_integer_
    gene_end   <- gc$end        %||% NA_integer_

    # ── GTEx v2: tissue-specific eQTL ─────────────────────────────────────
    # GTEx v8 uses versioned GENCODE v26 IDs (e.g. ENSG00000162434.11).
    # The Ensembl GRCh37 cache stores unversioned IDs; gtexr::get_genes()
    # is used once per gene to resolve the correct versioned ID.
    gtex_tissues_found <- character(0)

    if (!is.na(ensembl_id)) {
      # Resolve versioned GTEx GENCODE ID (cached under "gtex_gencode_id")
      gtex_gencode_id <- tryCatch({
        gid_entry <- failmode_cache_read("gtex_gencode_id", g, diag_cfg)
        if (!is.null(gid_entry)) {
          gid_entry$gencode_id
        } else {
          gid_tbl <- gtexr::get_genes(
            geneId         = toupper(g),
            gencodeVersion = "v26",
            genomeBuild    = "GRCh38/hg38"
          )
          gid <- if (!is.null(gid_tbl) && nrow(gid_tbl) > 0)
            gid_tbl$gencodeId[1] else ensembl_id
          failmode_cache_write(list(gencode_id = gid), "gtex_gencode_id", g, diag_cfg)
          gid
        }
      }, error = function(e) {
        failmode_log(log_file,
          "[fail_mode] gtex_gencode_id {g} lookup failed: {conditionMessage(e)} — using unversioned ID")
        ensembl_id
      })

      gtex_entry <- failmode_cache_read("gtex", g, diag_cfg)
      if (!is.null(gtex_entry)) {
        failmode_api_log_append(api_log, g, "gtex",
          "https://gtexportal.org/api/v2/association/singleTissueEqtl",
          cache_hit = TRUE, success = TRUE)
      } else {
        tissue_hits <- vapply(gtex_tissues, function(tissue) {
          url <- sprintf(
            "https://gtexportal.org/api/v2/association/singleTissueEqtl?gencodeId=%s&tissueSiteDetailId=%s&datasetId=gtex_v8&page=0&itemsPerPage=1",
            gtex_gencode_id, tissue)
          resp <- failmode_timed_get(url, rest_cfg, "gtex", g, api_log)
          hits <- resp$data %||% list()
          # Any hit returned means there is a significant eQTL in this tissue
          length(hits) > 0
        }, logical(1))
        gtex_entry <- list(tissue_hits = gtex_tissues[tissue_hits],
                           gtex_gencode_id = gtex_gencode_id)
        failmode_cache_write(gtex_entry, "gtex", g, diag_cfg)
      }
      gtex_tissues_found <- gtex_entry$tissue_hits %||% character(0)
    }

    gtex_blood_eqtl <- "Whole_Blood" %in% gtex_tissues_found
    gtex_eqtl_count <- length(gtex_tissues_found)

    # ── eQTL Catalogue v2: eQTL (cell-type) + sQTL ───────────────────────
    # Endpoint:  /v2/datasets/{QTD_ID}/associations?gene_id={ENSEMBL}&nlog10p=5&size=10
    # 400 "No results" → failmode_timed_get returns NULL → treated as 0 hits (correct)
    # All returned hits already pass nlog10p≥5 (p<1e-5); no re-filtering needed.
    eqtl_cell_hits <- character(0)
    sqtl_found     <- FALSE

    if (!is.na(ensembl_id)) {
      eqtl_entry <- failmode_cache_read("eqtl_catalogue", g, diag_cfg)
      if (!is.null(eqtl_entry)) {
        failmode_api_log_append(api_log, g, "eqtl_catalogue",
          "https://www.ebi.ac.uk/eqtl/api/v2/datasets/{QTD}/associations",
          cache_hit = TRUE, success = TRUE)
      } else {
        # ge eQTL: check immune cell type datasets
        for (ds in eqtl_ge_datasets) {
          url_e <- sprintf(
            "https://www.ebi.ac.uk/eqtl/api/v2/datasets/%s/associations?gene_id=%s&nlog10p=5&size=10",
            ds, ensembl_id)
          resp_e <- failmode_timed_get(url_e, eqtl_cfg, "eqtl_catalogue", g, api_log)
          # v2 response is a raw JSON array (unnamed list); NULL = HTTP 400 "No results"
          hits_e <- if (is.null(resp_e)) list()
                    else if (is.list(resp_e) && is.null(names(resp_e))) resp_e
                    else resp_e$results %||% list()
          if (length(hits_e) > 0) eqtl_cell_hits <- union(eqtl_cell_hits, ds)
        }

        # sQTL: leafcutter / txrev immune datasets
        for (ds in eqtl_sqtl_datasets) {
          url_s <- sprintf(
            "https://www.ebi.ac.uk/eqtl/api/v2/datasets/%s/associations?gene_id=%s&nlog10p=5&size=10",
            ds, ensembl_id)
          resp_s <- failmode_timed_get(url_s, eqtl_cfg, "eqtl_catalogue", g, api_log)
          hits_s <- if (is.null(resp_s)) list()
                    else if (is.list(resp_s) && is.null(names(resp_s))) resp_s
                    else resp_s$results %||% list()
          if (length(hits_s) > 0) sqtl_found <- TRUE
        }

        eqtl_entry <- list(cell_hits = eqtl_cell_hits, sqtl_found = sqtl_found)
        failmode_cache_write(eqtl_entry, "eqtl_catalogue", g, diag_cfg)
      }
      eqtl_cell_hits <- eqtl_entry$cell_hits  %||% character(0)
      sqtl_found     <- isTRUE(eqtl_entry$sqtl_found)
    }

    aries_entry <- failmode_cache_read("aries_mqtl", g, diag_cfg)
    if (!is.null(aries_entry)) {
      mqtl_result <- aries_entry
      failmode_api_log_append(api_log, g, "mqtl",
        paste0("cache/aries_mqtl/", toupper(g)),
        cache_hit = TRUE, success = TRUE,
        note = paste0("mqtl_status=", mqtl_result$mqtl_status,
                      " n_hits=", mqtl_result$n_hits %||% NA))
    } else if (!is.null(aries_data)) {
      mqtl_result <- fetch_aries_mqtl(g, aries_data, pval_threshold = aries_pval_thr)
      failmode_cache_write(mqtl_result, "aries_mqtl", g, diag_cfg)
      failmode_api_log_append(api_log, g, "mqtl",
        paste0("aries_mrinstruments/", toupper(g)),
        cache_hit = FALSE, success = TRUE, http_status = NA_integer_,
        note = paste0("mqtl_status=", mqtl_result$mqtl_status,
                      " n_hits=", mqtl_result$n_hits,
                      " top_cpg=", mqtl_result$top_cpg %||% "none"))
    } else {
      mqtl_result <- list(
        mqtl_status = "data_unavailable", n_hits = NA_integer_,
        top_snp = NA_character_, top_cpg = NA_character_, top_p = NA_real_,
        source = "aries_load_failed",
        note = "MRInstruments::aries_mqtl could not be loaded"
      )
      failmode_api_log_append(api_log, g, "mqtl",
        "aries_load_failed",
        cache_hit = FALSE, success = FALSE,
        note = mqtl_result$note)
    }
    godmc_has_cis_mqtl <- identical(mqtl_result$mqtl_status, "tested_positive")
    mqtl_status        <- mqtl_result$mqtl_status

    # ── Chromatin state (Roadmap/ENCODE ChromHMM) ─────────────────────────
    chr_blood    <- NA_character_
    chr_tcell    <- NA_character_
    chr_bcell    <- NA_character_
    chr_monocyte <- NA_character_

    if (chr_enabled) {
      chr_entry <- failmode_cache_read("chromatin", g, diag_cfg)
      if (!is.null(chr_entry)) {
        chr_blood    <- chr_entry$blood
        chr_tcell    <- chr_entry$tcell
        chr_bcell    <- chr_entry$bcell
        chr_monocyte <- chr_entry$monocyte
      } else if (!is.null(chromatin_beds) && !is.na(gene_chr) && !is.na(gene_start)) {
        chr_blood    <- parse_chromhmm_label(
          get_chromatin_state_at_tss(gene_chr, gene_start, chromatin_beds$blood,    tss_window))
        chr_tcell    <- parse_chromhmm_label(
          get_chromatin_state_at_tss(gene_chr, gene_start, chromatin_beds$tcell,    tss_window))
        chr_bcell    <- parse_chromhmm_label(
          get_chromatin_state_at_tss(gene_chr, gene_start, chromatin_beds$bcell,    tss_window))
        chr_monocyte <- parse_chromhmm_label(
          get_chromatin_state_at_tss(gene_chr, gene_start, chromatin_beds$monocyte, tss_window))
        failmode_cache_write(
          list(blood = chr_blood, tcell = chr_tcell,
               bcell = chr_bcell, monocyte = chr_monocyte),
          "chromatin", g, diag_cfg
        )
      }
    }

    # ── QTLbase2: local pre-queried xQTL ──────────────────────────────────
    q2g <- qtlbase2_data[[g]] %||% list()

    list(gene = g,
         gtex_eqtl_tissues         = gtex_tissues_found,
         gtex_eqtl_count           = gtex_eqtl_count,
         gtex_blood_eqtl           = gtex_blood_eqtl,
         eqtl_cat_cell_types       = eqtl_cell_hits,
         eqtl_cat_sqtl_hits        = sqtl_found,
         godmc_has_cis_mqtl        = godmc_has_cis_mqtl,   # NA = untested
         mqtl_status               = mqtl_status,          # "data_unavailable" | "tested_positive" | "tested_negative" | "api_error" | "no_instrument"
         qtlbase2_has_mqtl         = isTRUE(q2g$has_mqtl),
         qtlbase2_has_hqtl         = isTRUE(q2g$has_hqtl),
         qtlbase2_has_caqtl        = isTRUE(q2g$has_caqtl),
         qtlbase2_has_sqtl         = isTRUE(q2g$has_sqtl),
         qtlbase2_xqtl_types       = q2g$xqtl_types %||% character(0),
         chromatin_state_blood     = chr_blood,
         chromatin_state_tcell     = chr_tcell,
         chromatin_state_bcell     = chr_bcell,
         chromatin_state_monocyte  = chr_monocyte,
         chromatin_active_marks    = list(active_mks),
         chromatin_repressed_marks = list(repressed_mks))
  })

  names(result) <- target_genes
  n_ok <- sum(vapply(result, function(x) x$gtex_eqtl_count > 0 || isTRUE(x$eqtl_cat_sqtl_hits), logical(1)))
  failmode_log(log_file,
    "[fail_mode] fetch_multiomics_qtl done genes_with_any_qtl={n_ok}/{length(target_genes)}")
  result
}


# ─────────────────────────────────────────────────────────────────────────────
#
# ─────────────────────────────────────────────────────────────────────────────

run_cat2_structural_variant <- function(internal_cache, genomic_cache, diag_cfg,
                                        log_file = NULL) {
  cfg2 <- diag_cfg$cat2_structural_variant %||% list()
  if (!isTRUE(cfg2$enabled %||% TRUE)) {
    failmode_log(log_file, "[fail_mode] cat2 disabled, skipping")
    return(data.frame())
  }

  target_genes <- identify_diagnostic_targets(internal_cache, diag_cfg)
  if (length(target_genes) == 0) return(data.frame())

  failmode_log(log_file, "[fail_mode] cat2 start n={length(target_genes)}")

  class_df <- internal_cache$classification

  known_clusters   <- diag_cfg$fetch_genomic_context$paralog_query$known_clusters %||% list()
  fcgr_genes       <- character(0)
  for (cl in known_clusters) {
    if (!is.null(cl$chromosome) && cl$chromosome == "1") {
      fcgr_genes <- union(fcgr_genes, tolower(cl$genes %||% character(0)))
    }
  }

  # Scoring weights + thresholds
  sw <- cfg2$complex_locus_scoring$weights %||%
    list(MHC_CORE = 3, MHC_EXTENDED = 2, CAT2_STRUCTURAL_CLUSTER = 2,
         STRUCTURAL_CLUSTER = 2, HIGH_SV_BURDEN = 1, MODERATE_SV_BURDEN = 0.5,
         LOW_MAPPABILITY = 1)
  st <- cfg2$complex_locus_scoring$confidence_thresholds %||%
    list(high = 3, moderate = 2, low = 1)
  min_par       <- cfg2$detection_rules$structural_cluster$parameters$min_paralogs_for_cluster %||% 2L
  high_sv       <- cfg2$detection_rules$high_sv_burden$parameters$high_sv_burden_threshold     %||% 5L
  low_map       <- cfg2$detection_rules$low_mappability$parameters$low_mappability_threshold    %||% 0.5
  high_id_thr   <- 80.0   # CAT2_HIGH_PARALOG_HOMOLOGY: ≥2 paralogs AND identity ≥80%
  mod_id_thr    <- 60.0   # CAT2_MODERATE_PARALOG:      ≥1 paralog  AND identity ≥60%

  rows <- lapply(target_genes, function(g) {
    class_row  <- if (!is.null(class_df) && "gene" %in% names(class_df))
      class_df[tolower(class_df$gene) == g, , drop = FALSE] else data.frame()
    gene_class <- failmode_get_chr(class_row, "class")

    gc                   <- genomic_cache[[g]] %||% list()
    n_paralogs_nearby    <- gc$n_paralogs_nearby    %||% NA_integer_
    n_paralogs_all       <- gc$n_paralogs_all       %||% NA_integer_
    paralogs_nearby_vec  <- gc$paralogs_nearby      %||% character(0)
    max_par_identity     <- gc$max_paralog_identity %||% NA_real_
    sv_count             <- gc$sv_count             %||% NA_real_
    map_mean             <- gc$mappability_mean     %||% NA_real_

    # Use total paralog count (all) if available, fall back to nearby
    n_paralogs <- if (!is.na(n_paralogs_all)) n_paralogs_all else n_paralogs_nearby

    # target_annotations fallback
    locus_flags_raw  <- diag_cfg$target_annotations[[g]]$known_locus_flags
    locus_flags      <- if (!is.null(locus_flags_raw)) tolower(locus_flags_raw) else character(0)
    has_struct_annot <- "structural_cluster" %in% locus_flags
    cat1_flag_pre    <- if ("mhc_core" %in% locus_flags) "MHC_CORE"
                        else if ("mhc_extended" %in% locus_flags) "MHC_EXTENDED"
                        else NA_character_

    # Known structural cluster check (YAML annotation OR FCGR gene)
    in_known_cluster <- has_struct_annot || g %in% fcgr_genes
    cluster_name     <- if (in_known_cluster) {
      cl_match <- Filter(function(cl) g %in% tolower(cl$genes %||% character(0)), known_clusters)
      if (length(cl_match) > 0) names(cl_match)[1] %||% "unknown_cluster" else "FCGR_cluster"
    } else NA_character_

    # Priority: CAT2_STRUCTURAL_CLUSTER > CAT2_HIGH_PARALOG_HOMOLOGY > CAT2_MODERATE_PARALOG
    cat2_flag <- if (in_known_cluster) {
      "CAT2_STRUCTURAL_CLUSTER"
    } else if (!is.na(n_paralogs) && n_paralogs >= 2L &&
               !is.na(max_par_identity) && max_par_identity >= high_id_thr) {
      "CAT2_HIGH_PARALOG_HOMOLOGY"
    } else if (!is.na(n_paralogs) && n_paralogs >= 1L &&
               !is.na(max_par_identity) && max_par_identity >= mod_id_thr) {
      "CAT2_MODERATE_PARALOG"
    } else {
      NA_character_
    }

    cat2_confidence <- if (identical(cat2_flag, "CAT2_STRUCTURAL_CLUSTER"))    "High"
    else if (identical(cat2_flag, "CAT2_HIGH_PARALOG_HOMOLOGY"))               "High"
    else if (identical(cat2_flag, "CAT2_MODERATE_PARALOG"))                    "Moderate"
    else                                                                        NA_character_

    # Supplementary SV/mappability info (no longer primary flags, kept for reference)
    sv_flag  <- if (!is.na(sv_count)) {
      if (sv_count >= high_sv) "HIGH_SV_BURDEN"
      else if (sv_count >= 3L) "MODERATE_SV_BURDEN"
      else NA_character_
    } else NA_character_
    map_flag <- if (!is.na(map_mean) && map_mean < low_map) "LOW_MAPPABILITY" else NA_character_

    # Composite locus complexity score (cat1 + cat2 flags, SV/map as supplementary)
    active_flags <- unique(c(cat1_flag_pre, cat2_flag, sv_flag, map_flag))
    active_flags <- active_flags[!is.na(active_flags)]
    locus_score  <- sum(vapply(active_flags, function(f) as.numeric(sw[[f]] %||% 0), numeric(1)))
    locus_conf   <- if (locus_score >= as.numeric(st$high))     "High"
    else if (locus_score >= as.numeric(st$moderate))            "Moderate"
    else if (locus_score >= as.numeric(st$low))                 "Low"
    else                                                        NA_character_

    alt_strategy <- lookup_alternative_strategy(cat2_flag, diag_cfg)

    data.frame(
      gene                        = g,
      class                       = gene_class,
      cat2_flag                   = cat2_flag,
      cat2_confidence             = cat2_confidence,
      cat2_locus_score            = locus_score,
      cat2_locus_score_confidence = locus_conf,
      paralog_count               = n_paralogs,
      max_paralog_identity        = max_par_identity,
      in_known_cluster            = in_known_cluster,
      cluster_name                = cluster_name,
      n_paralogs_nearby           = n_paralogs_nearby,
      paralogs_nearby             = paste(paralogs_nearby_vec, collapse = ";"),
      sv_count                    = sv_count,
      mappability_mean            = map_mean,
      ensembl_id                  = gc$ensembl_id %||% NA_character_,
      chr                         = gc$chr        %||% NA_character_,
      recommended_strategy        = alt_strategy,
      annotation_timestamp        = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
      stringsAsFactors            = FALSE
    )
  })

  result <- dplyr::bind_rows(rows)
  n_flagged <- sum(!is.na(result$cat2_flag), na.rm = TRUE)
  failmode_log(log_file,
    "[fail_mode] cat2 done n_genes={nrow(result)} n_flagged={n_flagged}")
  result
}


# ─────────────────────────────────────────────────────────────────────────────
#
#
#   pQTL  — internal_cache$pqtl_preferred (pipeline_internal)
#   mQTL  — godmc_has_cis_mqtl OR qtlbase2_has_mqtl/hqtl/caqtl
#   eQTL  — gtex_blood_eqtl
# ─────────────────────────────────────────────────────────────────────────────

run_cat4_tissue_cell_specificity <- function(internal_cache, qtl_cache, diag_cfg,
                                             log_file = NULL) {
  cfg4 <- diag_cfg$cat4_tissue_cell_specificity %||% list()
  if (!isTRUE(cfg4$enabled %||% TRUE)) {
    failmode_log(log_file, "[fail_mode] cat4 disabled, skipping")
    return(data.frame())
  }

  target_genes <- identify_diagnostic_targets(internal_cache, diag_cfg)
  if (length(target_genes) == 0) return(data.frame())

  failmode_log(log_file, "[fail_mode] cat4 start n={length(target_genes)}")

  class_df  <- internal_cache$classification
  pqtl_pref <- internal_cache$pqtl_preferred

  rows <- lapply(target_genes, function(g) {
    class_row  <- if (!is.null(class_df) && "gene" %in% names(class_df))
      class_df[tolower(class_df$gene) == g, , drop = FALSE] else data.frame()
    gene_class <- failmode_get_chr(class_row, "class")

    qc <- qtl_cache[[g]] %||% list()

    sqtl_detected <- isTRUE(qc$eqtl_cat_sqtl_hits) || isTRUE(qc$qtlbase2_has_sqtl)

    # ── mQTL available (GoDMC OR QTLbase2) ────────────────────────────────
    mqtl_status_g  <- qc$mqtl_status %||% "data_unavailable"
    qtlbase2_mqtl  <- isTRUE(qc$qtlbase2_has_mqtl) ||
                      isTRUE(qc$qtlbase2_has_hqtl)  ||
                      isTRUE(qc$qtlbase2_has_caqtl)
    mqtl_available <- if (qtlbase2_mqtl || mqtl_status_g == "tested_positive") {
      TRUE
    } else if (mqtl_status_g == "tested_negative") {
      FALSE
    } else {
    }

    # ── pQTL available (pipeline_internal pqtl_preferred) ─────────────────
    pqtl_available <- !is.null(pqtl_pref) &&
                      is.data.frame(pqtl_pref) &&
                      nrow(pqtl_pref) > 0 &&
                      "gene" %in% names(pqtl_pref) &&
                      any(tolower(pqtl_pref$gene) == g)

    # ── eQTL significant (GTEx blood eQTL) ────────────────────────────────
    eqtl_significant <- isTRUE(qc$gtex_blood_eqtl)

    # ── Cat4 flag: CAT4_MULTIOMICS_CONVERGENCE > CAT4_EPIGENETIC_CONTEXT ──
    cat4_flag <- if (sqtl_detected && pqtl_available) {
      "CAT4_MULTIOMICS_CONVERGENCE"
    } else if (mqtl_available && eqtl_significant) {
      "CAT4_EPIGENETIC_CONTEXT"
    } else {
      NA_character_
    }

    cat4_confidence <- if (identical(cat4_flag, "CAT4_MULTIOMICS_CONVERGENCE")) "High"
    else if (identical(cat4_flag, "CAT4_EPIGENETIC_CONTEXT"))                   "Moderate"
    else                                                                         NA_character_

    alt_strategy <- lookup_alternative_strategy(cat4_flag, diag_cfg)

    data.frame(
      gene                 = g,
      class                = gene_class,
      cat4_flag            = cat4_flag,
      cat4_confidence      = cat4_confidence,
      sqtl_detected        = sqtl_detected,
      mqtl_available       = mqtl_available,  # TRUE/FALSE/NA (NA = data_unavailable)
      mqtl_status          = mqtl_status_g,
      pqtl_available       = pqtl_available,
      eqtl_significant     = eqtl_significant,
      recommended_strategy = alt_strategy,
      annotation_timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
      stringsAsFactors     = FALSE
    )
  })

  result <- dplyr::bind_rows(rows)
  n_flagged <- sum(!is.na(result$cat4_flag), na.rm = TRUE)
  failmode_log(log_file,
    "[fail_mode] cat4 done n_genes={nrow(result)} n_flagged={n_flagged}")
  result
}
