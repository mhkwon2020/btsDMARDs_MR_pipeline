# R/07_mr_estimation.R
# MR estimation utilities
# Version: 1.0
# Author: Mihye Kwon
suppressPackageStartupMessages({
  library(dplyr)
  library(TwoSampleMR)
})

get_egger_minimum_nsnp <- function(protocol_cfg = NULL) {
  if (!is.null(protocol_cfg) &&
      !is.null(protocol_cfg$mr_estimation$method_rules$mr_egger$minimum_nsnp)) {
    return(as.integer(protocol_cfg$mr_estimation$method_rules$mr_egger$minimum_nsnp))
  }

  if (!is.null(protocol_cfg) &&
      !is.null(protocol_cfg$instrument_selection$minimum_instrument_rules$require_three_snp_for_egger)) {
    return(if (isTRUE(protocol_cfg$instrument_selection$minimum_instrument_rules$require_three_snp_for_egger)) 3L else 2L)
  }

  3L
}

safe_calc_i2_gx <- function(usable_harm) {
  if (is.null(usable_harm) || nrow(usable_harm) < 2) return(NA_real_)
  if (!all(c("beta.exposure", "se.exposure") %in% names(usable_harm))) return(NA_real_)

  tryCatch(
    {
      y <- abs(as.numeric(usable_harm$beta.exposure))
      s <- as.numeric(usable_harm$se.exposure)
      if (all(is.finite(y)) && all(is.finite(s))) {
        TwoSampleMR::Isq(y, s)
      } else {
        NA_real_
      }
    },
    error = function(e) NA_real_
  )
}

classify_egger_reliability <- function(i2_gx, nsnp_usable, egger_min) {
  if (is.na(i2_gx)) return("unknown")
  if (nsnp_usable < egger_min) return("not_applicable")
  if (nsnp_usable < 10) return("low_information_nsnp_lt_10")
  if (i2_gx < 0.90) return("possible_nome_violation")
  "acceptable"
}

select_mr_methods_by_nsnp <- function(nsnp_usable, protocol_cfg = NULL) {
  if (nsnp_usable <= 0) return(character(0))
  if (nsnp_usable == 1) return("mr_wald_ratio")
  if (nsnp_usable == 2) return("mr_ivw")

  egger_min <- get_egger_minimum_nsnp(protocol_cfg)

  methods <- c("mr_ivw", "mr_weighted_median")
  if (nsnp_usable >= egger_min) {
    methods <- c(methods, "mr_egger_regression")
  }
  methods
}

summarise_f_statistics <- function(usable_harm) {
  if (is.null(usable_harm) || nrow(usable_harm) == 0) {
    return(data.frame(
      f_median = NA_real_,
      f_min = NA_real_,
      f_max = NA_real_,
      n_f_gt_10 = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }

  if (!all(c("beta.exposure", "se.exposure") %in% names(usable_harm))) {
    return(data.frame(
      f_median = NA_real_,
      f_min = NA_real_,
      f_max = NA_real_,
      n_f_gt_10 = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }

  fvec <- (usable_harm$beta.exposure / usable_harm$se.exposure)^2
  fvec <- fvec[is.finite(fvec)]

  if (length(fvec) == 0) {
    return(data.frame(
      f_median = NA_real_,
      f_min = NA_real_,
      f_max = NA_real_,
      n_f_gt_10 = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    f_median = stats::median(fvec, na.rm = TRUE),
    f_min = min(fvec, na.rm = TRUE),
    f_max = max(fvec, na.rm = TRUE),
    n_f_gt_10 = sum(fvec > 10, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

run_mr_suite <- function(harm, outcome_obj, run_cfg, protocol_cfg = NULL,
                         log_file = NULL, err_file = NULL) {
  if (is.null(harm) || nrow(harm) == 0) {
    return(list(
      mr = NULL,
      heterogeneity = NULL,
      pleiotropy = NULL,
      loo = NULL,
      f_summary = data.frame()
    ))
  }

  usable <- harm[harm$mr_keep, , drop = FALSE]
  nsnp_usable <- length(unique(usable$SNP))

  methods <- select_mr_methods_by_nsnp(nsnp_usable, protocol_cfg = protocol_cfg)
  if (length(methods) == 0) {
    return(list(
      mr = NULL,
      heterogeneity = NULL,
      pleiotropy = NULL,
      loo = NULL,
      f_summary = summarise_f_statistics(usable)
    ))
  }

  f_summary <- summarise_f_statistics(usable)

  mrres <- tryCatch(
    TwoSampleMR::mr(usable, method_list = methods),
    error = function(e) {
      if (!is.null(err_file)) {
        log_message(err_file, "MR FAIL | outcome={outcome_obj$name} gene={unique(usable$exposure)} msg={e$message}")
      }
      NULL
    }
  )

  if (!is.null(mrres) && nrow(mrres) > 0) {
    mrres <- mrres %>%
      mutate(
        gene = unique(usable$exposure),
        outcome_name = outcome_obj$name,
        outcome_id = outcome_obj$id,
        role = outcome_obj$role,
        ancestry = outcome_obj$ancestry,
        phenotype = outcome_obj$phenotype %||% NA_character_,
        nsnp_usable = nsnp_usable,
        OR = exp(b),
        OR_lci95 = exp(b - 1.96 * se),
        OR_uci95 = exp(b + 1.96 * se)
      ) %>%
      bind_cols(f_summary[rep(1, nrow(.)), , drop = FALSE])
  }

  het <- NULL
  if (isTRUE(run_cfg$sensitivity_to_run$heterogeneity) && nsnp_usable >= 2) {
    het <- tryCatch(
      TwoSampleMR::mr_heterogeneity(usable),
      error = function(e) NULL
    )
    if (!is.null(het) && nrow(het) > 0) {
      het$gene <- unique(usable$exposure)
      het$outcome_name <- outcome_obj$name
      het$outcome_id <- outcome_obj$id
      het$role <- outcome_obj$role
      het$ancestry <- outcome_obj$ancestry
      if ("Q" %in% names(het)) het$Q_statistic <- het$Q
      if ("Q_df" %in% names(het)) het$Q_df_statistic <- het$Q_df
      if ("Q_pval" %in% names(het)) het$Q_pvalue <- het$Q_pval
    }
  }

  ple <- NULL
  egger_min <- get_egger_minimum_nsnp(protocol_cfg)

  if (isTRUE(run_cfg$sensitivity_to_run$pleiotropy) && nsnp_usable >= egger_min) {
    ple <- tryCatch(
      TwoSampleMR::mr_pleiotropy_test(usable),
      error = function(e) NULL
    )
    if (!is.null(ple) && nrow(ple) > 0) {
      i2_gx <- if (isTRUE(run_cfg$sensitivity_to_run$report_i2_gx_with_egger %||% TRUE)) {
        safe_calc_i2_gx(usable)
      } else {
        NA_real_
      }
      ple$gene <- unique(usable$exposure)
      ple$outcome_name <- outcome_obj$name
      ple$outcome_id <- outcome_obj$id
      ple$role <- outcome_obj$role
      ple$ancestry <- outcome_obj$ancestry
      ple$egger_minimum_nsnp <- egger_min
      ple$i2_gx <- i2_gx
      ple$egger_reliability <- classify_egger_reliability(i2_gx, nsnp_usable, egger_min)
      ple$nsnp_usable <- nsnp_usable
    }
  }

  loo <- NULL
  if (isTRUE(run_cfg$sensitivity_to_run$leave_one_out) && nsnp_usable >= 3) {
    loo <- tryCatch(
      TwoSampleMR::mr_leaveoneout(usable),
      error = function(e) NULL
    )
    if (!is.null(loo) && nrow(loo) > 0) {
      loo$gene <- unique(usable$exposure)
      loo$outcome_name <- outcome_obj$name
    }
  }

  list(
    mr = mrres,
    heterogeneity = het,
    pleiotropy = ple,
    loo = loo,
    f_summary = f_summary
  )
}
