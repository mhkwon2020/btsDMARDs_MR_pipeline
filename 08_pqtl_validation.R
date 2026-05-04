# R/08_pqtl_validation.R
# pQTL validation utilities
# Version: 1.0
# Author: Mihye Kwon
suppressPackageStartupMessages({
  library(dplyr)
})

get_pqtl_outcome_scope_key <- function(outcome_obj) {
  if (identical(outcome_obj$role %||% NA_character_, "primary")) return("primary")
  if (identical(outcome_obj$role %||% NA_character_, "replication")) return("replication")
  if (identical(outcome_obj$name %||% NA_character_, "RA_EUR_SEROPOS")) return("seropos")
  if (identical(outcome_obj$name %||% NA_character_, "RA_EUR_SERONEG")) return("seroneg")
  if (identical(outcome_obj$name %||% NA_character_, "RA_EAS_BBJ")) return("eas")
  tolower(outcome_obj$role %||% outcome_obj$name %||% "unknown")
}

outcome_allowed_for_pqtl_source <- function(outcome_obj, src) {
  allowed <- src$outcome_scope %||% character(0)
  if (length(allowed) == 0) return(TRUE)
  scope_key <- get_pqtl_outcome_scope_key(outcome_obj)
  scope_key %in% tolower(as.character(allowed))
}

build_preferred_results_by_source <- function(all_mr_results, run_cfg, source_col = "pqtl_source") {
  if (is.null(all_mr_results) || nrow(all_mr_results) == 0) return(data.frame())

  method_order <- preferred_method_order(run_cfg)
  split_key <- interaction(
    all_mr_results$gene,
    all_mr_results$outcome_name,
    all_mr_results[[source_col]],
    drop = TRUE
  )

  split_df <- split(all_mr_results, split_key, drop = TRUE)
  preferred <- lapply(split_df, function(df) pick_preferred_result(df, method_order))
  bind_rows(preferred)
}

run_pqtl_validation_suite <- function(pqtl_sources,
                                      outcomes,
                                      genes,
                                      protocol_cfg,
                                      run_cfg,
                                      target_meta = NULL,
                                      log_file = NULL,
                                      err_file = NULL) {
  all_attr <- list()
  all_mr <- list()
  all_het <- list()
  all_ple <- list()
  all_loo <- list()

  if (length(pqtl_sources) == 0) {
    return(list(
      attrition = data.frame(),
      mr = data.frame(),
      heterogeneity = data.frame(),
      pleiotropy = data.frame(),
      loo = data.frame(),
      preferred = data.frame()
    ))
  }

  for (src in pqtl_sources) {
    exp_raw <- src$data
    if (is.null(exp_raw) || nrow(exp_raw) == 0) next

    for (o in outcomes) {
      if (!outcome_allowed_for_pqtl_source(o, src)) {
        if (!is.null(log_file)) {
          log_message(log_file, "PQTLSOURCE skip source={src$source_name} outcome={o$name} reason=outcome_scope_filtered")
        }
        next
      }

      pop <- toupper(o$ancestry %||% "EUR")
      if (!is.null(log_file)) {
        log_message(log_file, "PQTLSOURCE start source={src$source_name} outcome={o$name} ancestry={pop}")
      }

      for (g in genes) {
        inst <- build_instrument_for_gene(
          exp_raw = exp_raw,
          gene = g,
          ancestry = pop,
          protocol_cfg = protocol_cfg,
          target_meta = target_meta,
          log_file = log_file,
          err_file = err_file
        )

        harm_res <- extract_and_harmonise(
          exp_dat = inst$exp_dat,
          outcome_obj = o,
          protocol_cfg = protocol_cfg,
          log_file = log_file,
          err_file = err_file
        )

        attr_row <- inst$attrition
        attr_row$post_outcome_extraction_count <- harm_res$attrition$post_outcome_extraction_count
        attr_row$post_harmonisation_row_count <- harm_res$attrition$post_harmonisation_row_count
        attr_row$post_harmonisation_usable_count <- harm_res$attrition$post_harmonisation_usable_count
        attr_row$palindromic_loss_count <- harm_res$attrition$palindromic_loss_count

        if (is.na(attr_row$non_estimable_reason) && !is.na(harm_res$attrition$non_estimable_reason)) {
          attr_row$non_estimable_reason <- harm_res$attrition$non_estimable_reason
        }

        attr_row$outcome_name <- o$name
        attr_row$outcome_id <- o$id
        attr_row$role <- o$role
        attr_row$phenotype <- o$phenotype %||% NA_character_
        attr_row$pqtl_source <- src$source_name
        attr_row$pqtl_source_id <- src$source_id
        attr_row$proxy_layer <- "pQTL"
        attr_row$estimable <- !is.null(harm_res$harm) &&
          nrow(harm_res$harm) > 0 &&
          sum(harm_res$harm$mr_keep, na.rm = TRUE) > 0

        if (is.na(attr_row$non_estimable_reason) && !isTRUE(attr_row$estimable)) {
          attr_row$non_estimable_reason <- "insufficient_variant_count_for_estimator"
        }

        all_attr[[length(all_attr) + 1]] <- attr_row

        if (isTRUE(attr_row$estimable)) {
          mr_suite <- run_mr_suite(
            harm = harm_res$harm,
            outcome_obj = o,
            run_cfg = run_cfg,
            protocol_cfg = protocol_cfg,
            log_file = log_file,
            err_file = err_file
          )

          if (!is.null(mr_suite$mr) && nrow(mr_suite$mr) > 0) {
            mr_piece <- mr_suite$mr
            mr_piece$pqtl_source <- src$source_name
            mr_piece$pqtl_source_id <- src$source_id
            mr_piece$proxy_layer <- "pQTL"
            all_mr[[length(all_mr) + 1]] <- mr_piece
          }
          if (!is.null(mr_suite$heterogeneity) && nrow(mr_suite$heterogeneity) > 0) {
            het_piece <- mr_suite$heterogeneity
            het_piece$pqtl_source <- src$source_name
            het_piece$pqtl_source_id <- src$source_id
            het_piece$proxy_layer <- "pQTL"
            all_het[[length(all_het) + 1]] <- het_piece
          }
          if (!is.null(mr_suite$pleiotropy) && nrow(mr_suite$pleiotropy) > 0) {
            ple_piece <- mr_suite$pleiotropy
            ple_piece$pqtl_source <- src$source_name
            ple_piece$pqtl_source_id <- src$source_id
            ple_piece$proxy_layer <- "pQTL"
            all_ple[[length(all_ple) + 1]] <- ple_piece
          }
          if (!is.null(mr_suite$loo) && nrow(mr_suite$loo) > 0) {
            loo_piece <- mr_suite$loo
            loo_piece$pqtl_source <- src$source_name
            loo_piece$pqtl_source_id <- src$source_id
            loo_piece$proxy_layer <- "pQTL"
            all_loo[[length(all_loo) + 1]] <- loo_piece
          }
        }
      }
    }
  }

  attrition_df <- if (length(all_attr) > 0) bind_rows(all_attr) else data.frame()
  mr_df <- if (length(all_mr) > 0) bind_rows(all_mr) else data.frame()
  het_df <- if (length(all_het) > 0) bind_rows(all_het) else data.frame()
  ple_df <- if (length(all_ple) > 0) bind_rows(all_ple) else data.frame()
  loo_df <- if (length(all_loo) > 0) bind_rows(all_loo) else data.frame()
  preferred_df <- build_preferred_results_by_source(mr_df, run_cfg, source_col = "pqtl_source")

  list(
    attrition = attrition_df,
    mr = mr_df,
    heterogeneity = het_df,
    pleiotropy = ple_df,
    loo = loo_df,
    preferred = preferred_df
  )
}
