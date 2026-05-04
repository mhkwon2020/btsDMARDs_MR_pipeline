# Main pipeline runner
# Version: 1.0
# Author: Mihye Kwon
#
# Dependencies:
# - Requires R/00_utils.R for %||% and log_message
# - Requires R/07_mr_estimation.R for run_mr_suite


suppressPackageStartupMessages({
  library(fs)
  library(glue)
  library(yaml)
  library(dplyr)
})

# -----------------------------------------------------------------------------
# Resolve external dependencies
# -----------------------------------------------------------------------------

plink_bin <- Sys.getenv("PLINK_BIN")
ld_ref_eur <- Sys.getenv("LD_REF_EUR")
ld_ref_eas <- Sys.getenv("LD_REF_EAS")

if (plink_bin == "") {
  stop("PLINK_BIN is not set. Please define it in .Renviron or environment.")
}
if (!file.exists(plink_bin)) {
  stop("PLINK_BIN path does not exist: ", plink_bin)
}

if (ld_ref_eur != "" && !dir.exists(ld_ref_eur)) {
  warning("LD_REF_EUR path not found: ", ld_ref_eur)
}

if (ld_ref_eas != "" && !dir.exists(ld_ref_eas)) {
  warning("LD_REF_EAS path not found: ", ld_ref_eas)
}

# -----------------------------------------------------------------------------
# Create run-specific output directory
# -----------------------------------------------------------------------------

run_dir <- file.path("runs", format(Sys.time(), "%Y%m%d_%H%M%S"))
dir.create(run_dir, recursive = TRUE)

message("Output directory: ", run_dir)

# Optional: attach to config
options(run_dir = run_dir)

# Save reproducibility info
writeLines(capture.output(sessionInfo()),
           file.path(run_dir, "session_info.txt"))

writeLines(as.character(Sys.time()),
           file.path(run_dir, "run_timestamp.txt"))

# -----------------------------------------------------------------------------
# Create run-specific output directory
# -----------------------------------------------------------------------------
run_dir <- file.path("runs", format(Sys.time(), "%Y%m%d_%H%M%S"))

dir.create(run_dir, recursive = TRUE)

message("Output directory: ", run_dir)

# -----------------------------------------------------------------------------
# Core modules
# -----------------------------------------------------------------------------
source("R/00_utils.R")
source("R/01_protocol_checks.R")
source("R/02_input_readers.R")

# -----------------------------------------------------------------------------
# Instrument + harmonisation
# -----------------------------------------------------------------------------
source("R/03_instrument_builder.R")
source("R/04_outcome_harmonise.R")
source("R/05_sample_overlap_check.R")

# -----------------------------------------------------------------------------
# Data layers
# -----------------------------------------------------------------------------
source("R/06_pqtl_loader.R")

# -----------------------------------------------------------------------------
# MR estimation
# -----------------------------------------------------------------------------
source("R/07_mr_estimation.R")
source("R/08_pqtl_validation.R")

# -----------------------------------------------------------------------------
# Downstream inference
# -----------------------------------------------------------------------------
source("R/09_colocalisation.R")
source("R/10_classification.R")

# -----------------------------------------------------------------------------
# Reporting
# -----------------------------------------------------------------------------
source("R/11_reporting_tables.R")
source("R/12_reporting_figures.R")

# -----------------------------------------------------------------------------
# Triangulation
# -----------------------------------------------------------------------------
source("R/13_triangulation.R")
source("R/14_additional_figures.R")

# -----------------------------------------------------------------------------
# Failure mode diagnostics
# -----------------------------------------------------------------------------
source("R/15_failure_mode_diagnostic.R")
source("R/16_failure_mode_cat6.R")

# -----------------------------------------------------------------------------
# Main execution function
# -----------------------------------------------------------------------------

run_eqtl_primary_suite <- function(exp_raw,
                                   outcomes,
                                   genes,
                                   protocol_cfg,
                                   run_cfg,
                                   target_meta = NULL,
                                   log_file = NULL,
                                   err_file = NULL) {

  all_attr <- list()
  all_mr   <- list()
  all_het  <- list()
  all_ple  <- list()
  all_loo  <- list()

  for (o in outcomes) {

    pop <- toupper(o$ancestry %||% "EUR")
    log_message(log_file, "OUTCOME start name={o$name} id={o$id} ancestry={pop}")

    for (g in genes) {

      log_message(log_file, "GENE start gene={g} outcome={o$name}")

      inst <- build_instrument_for_gene(
        exp_raw      = exp_raw,
        gene         = g,
        ancestry     = pop,
        protocol_cfg = protocol_cfg,
        target_meta  = target_meta,
        log_file     = log_file,
        err_file     = err_file
      )

      harm_res <- extract_and_harmonise(
        exp_dat      = inst$exp_dat,
        outcome_obj  = o,
        protocol_cfg = protocol_cfg,
        log_file     = log_file,
        err_file     = err_file
      )

      attr_row <- inst$attrition
      attr_row$post_outcome_extraction_count      <- harm_res$attrition$post_outcome_extraction_count
      attr_row$post_harmonisation_row_count       <- harm_res$attrition$post_harmonisation_row_count
      attr_row$post_harmonisation_usable_count    <- harm_res$attrition$post_harmonisation_usable_count
      attr_row$palindromic_loss_count             <- harm_res$attrition$palindromic_loss_count

      if (is.na(attr_row$non_estimable_reason) &&
          !is.na(harm_res$attrition$non_estimable_reason)) {
        attr_row$non_estimable_reason <- harm_res$attrition$non_estimable_reason
      }

      attr_row$outcome_name <- o$name
      attr_row$outcome_id   <- o$id
      attr_row$role         <- o$role
      attr_row$phenotype    <- o$phenotype %||% NA_character_
      attr_row$proxy_layer  <- "eQTL"

      attr_row$estimable <- !is.null(harm_res$harm) &&
        nrow(harm_res$harm) > 0 &&
        sum(harm_res$harm$mr_keep, na.rm = TRUE) > 0

      if (is.na(attr_row$non_estimable_reason) && !isTRUE(attr_row$estimable)) {
        attr_row$non_estimable_reason <- "insufficient_variant_count_for_estimator"
      }

      all_attr[[length(all_attr) + 1]] <- attr_row

      if (isTRUE(attr_row$estimable)) {

        mr_suite <- run_mr_suite(
          harm         = harm_res$harm,
          outcome_obj  = o,
          run_cfg      = run_cfg,
          protocol_cfg = protocol_cfg,
          log_file     = log_file,
          err_file     = err_file
        )

        if (!is.null(mr_suite$mr) && nrow(mr_suite$mr) > 0) {
          mr_piece <- mr_suite$mr
          mr_piece$proxy_layer <- "eQTL"
          all_mr[[length(all_mr) + 1]] <- mr_piece
        }

        if (!is.null(mr_suite$heterogeneity) && nrow(mr_suite$heterogeneity) > 0) {
          het_piece <- mr_suite$heterogeneity
          het_piece$proxy_layer <- "eQTL"
          all_het[[length(all_het) + 1]] <- het_piece
        }

        if (!is.null(mr_suite$pleiotropy) && nrow(mr_suite$pleiotropy) > 0) {
          ple_piece <- mr_suite$pleiotropy
          ple_piece$proxy_layer <- "eQTL"
          all_ple[[length(all_ple) + 1]] <- ple_piece
        }

        if (!is.null(mr_suite$loo) && nrow(mr_suite$loo) > 0) {
          loo_piece <- mr_suite$loo
          loo_piece$proxy_layer <- "eQTL"
          all_loo[[length(all_loo) + 1]] <- loo_piece
        }
      }
    }
  }

  list(
    attrition     = dplyr::bind_rows(all_attr),
    mr            = dplyr::bind_rows(all_mr),
    heterogeneity = dplyr::bind_rows(all_het),
    pleiotropy    = dplyr::bind_rows(all_ple),
    loo           = dplyr::bind_rows(all_loo)
  )
}
