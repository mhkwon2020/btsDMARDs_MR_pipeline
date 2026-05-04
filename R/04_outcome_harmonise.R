# R/04_outcome_harmonise.R
# Outcome harmonisation utilities
# Version: 1.0
# Author: Mihye Kwon
suppressPackageStartupMessages({
  library(dplyr)
  library(TwoSampleMR)
})

extract_and_harmonise <- function(exp_dat, outcome_obj, protocol_cfg,
                                  log_file = NULL, err_file = NULL) {
  attr <- data.frame(
    post_outcome_extraction_count = NA_integer_,
    post_harmonisation_row_count = NA_integer_,
    post_harmonisation_usable_count = NA_integer_,
    palindromic_loss_count = NA_integer_,
    non_estimable_reason = NA_character_,
    stringsAsFactors = FALSE
  )

  if (is.null(exp_dat) || nrow(exp_dat) == 0) {
    attr$non_estimable_reason <- "no_variant_after_clumping"
    return(list(out_dat = NULL, harm = NULL, attrition = attr))
  }

  snps <- unique(exp_dat$SNP)

  out_dat <- tryCatch(
    TwoSampleMR::extract_outcome_data(snps = snps, outcomes = outcome_obj$id),
    error = function(e) {
      if (!is.null(err_file)) {
        log_message(err_file, "OUTCOME EXTRACT FAIL | outcome={outcome_obj$name} gene={unique(exp_dat$exposure)} msg={e$message}")
      }
      NULL
    }
  )

  attr$post_outcome_extraction_count <- if (is.null(out_dat)) 0L else length(unique(out_dat$SNP))

  if (is.null(out_dat) || nrow(out_dat) == 0) {
    attr$non_estimable_reason <- "no_variant_after_outcome_extraction"
    return(list(out_dat = NULL, harm = NULL, attrition = attr))
  }

  harm <- tryCatch(
    TwoSampleMR::harmonise_data(exp_dat, out_dat, action = protocol_cfg$harmonisation$action),
    error = function(e) {
      if (!is.null(err_file)) {
        log_message(err_file, "HARMONISE FAIL | outcome={outcome_obj$name} gene={unique(exp_dat$exposure)} msg={e$message}")
      }
      NULL
    }
  )

  attr$post_harmonisation_row_count <- if (is.null(harm)) 0L else nrow(harm)

  if (is.null(harm) || nrow(harm) == 0) {
    attr$post_harmonisation_usable_count <- 0L
    attr$non_estimable_reason <- "no_usable_variant_after_harmonisation"
    return(list(out_dat = out_dat, harm = NULL, attrition = attr))
  }

  usable_n <- sum(harm$mr_keep, na.rm = TRUE)
  attr$post_harmonisation_usable_count <- usable_n

  if ("palindromic" %in% names(harm)) {
    attr$palindromic_loss_count <- sum(harm$palindromic & !harm$mr_keep, na.rm = TRUE)
  } else {
    attr$palindromic_loss_count <- NA_integer_
  }

  if (usable_n == 0) {
    attr$non_estimable_reason <- "no_usable_variant_after_harmonisation"
  }

  list(out_dat = out_dat, harm = harm, attrition = attr)
}
