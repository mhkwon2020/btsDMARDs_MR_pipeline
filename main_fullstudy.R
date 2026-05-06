# Main pipeline runner
# Version: 1.0
# Author: Mihye Kwon
#
# USAGE:
#   Set working directory to the project root, then:
#     source("main_fullstudy.R")
#   Or from terminal:
#     Rscript main_fullstudy.R
#
# ENVIRONMENT VARIABLES (set in .Renviron before running):
#   PLINK_BIN       - path to PLINK 1.9 binary
#   LD_REF_EUR      - path prefix for 1000G EUR reference (.bed/.bim/.fam)
#   LD_REF_EAS      - path prefix for 1000G EAS reference (.bed/.bim/.fam)
#   INPUT_DATA_PATH - full path to eQTLGen exposure file (.txt)
#   PQTL_DATA_DIR   - directory containing pQTL data files
#   OPENGWAS_JWT    - JWT token for OpenGWAS API access

suppressPackageStartupMessages({
  library(fs)
  library(glue)
  library(yaml)
  library(dplyr)
})

# -----------------------------------------------------------------------------
# Working directory guard
# -----------------------------------------------------------------------------
if (!file.exists("configs/protocol_v1.yaml") || !file.exists("R/00_utils.R")) {
  stop(
    "Working directory must be the project root.\n",
    "Current directory: ", getwd(), "\n",
    "Please run setwd('path/to/btsDMARDs_MR_pipeline') before sourcing.",
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# Resolve external dependencies
# -----------------------------------------------------------------------------
plink_bin  <- Sys.getenv("PLINK_BIN")
ld_ref_eur <- Sys.getenv("LD_REF_EUR")
ld_ref_eas <- Sys.getenv("LD_REF_EAS")

if (plink_bin == "") stop("PLINK_BIN is not set. Please define it in .Renviron.")
if (!file.exists(plink_bin)) stop("PLINK_BIN path does not exist: ", plink_bin)
if (ld_ref_eur != "" && !file.exists(paste0(ld_ref_eur, ".bed"))) warning("LD_REF_EUR file not found: ", ld_ref_eur)
if (ld_ref_eas != "" && !file.exists(paste0(ld_ref_eas, ".bed"))) warning("LD_REF_EAS file not found: ", ld_ref_eas)

# -----------------------------------------------------------------------------
# Source all modules
# -----------------------------------------------------------------------------
source("R/00_utils.R")
source("R/01_protocol_checks.R")
source("R/02_input_readers.R")
source("R/03_instrument_builder.R")
source("R/04_outcome_harmonise.R")
source("R/05_sample_overlap_check.R")
source("R/06_pqtl_loader.R")
source("R/07_mr_estimation.R")
source("R/08_pqtl_validation.R")
source("R/09_colocalisation.R")
source("R/10_classification.R")
source("R/11_reporting_tables.R")
source("R/12_reporting_figures.R")
source("R/13_triangulation.R")
source("R/14_additional_figures.R")
source("R/15_failure_mode_diagnostic.R")
source("R/16_failure_mode_cat6.R")

# -----------------------------------------------------------------------------
# eQTL primary suite
# -----------------------------------------------------------------------------
run_eqtl_primary_suite <- function(exp_raw, outcomes, genes,
                                   protocol_cfg, run_cfg,
                                   target_meta = NULL,
                                   log_file = NULL, err_file = NULL) {
  all_attr <- list(); all_mr <- list(); all_het <- list()
  all_ple  <- list(); all_loo <- list()

  for (o in outcomes) {
    pop <- toupper(o$ancestry %||% "EUR")
    log_message(log_file, "OUTCOME start name={o$name} id={o$id} ancestry={pop}")

    for (g in genes) {
      log_message(log_file, "GENE start gene={g} outcome={o$name}")

      inst <- build_instrument_for_gene(
        exp_raw=exp_raw, gene=g, ancestry=pop,
        protocol_cfg=protocol_cfg, target_meta=target_meta,
        log_file=log_file, err_file=err_file)

      harm_res <- extract_and_harmonise(
        exp_dat=inst$exp_dat, outcome_obj=o,
        protocol_cfg=protocol_cfg, log_file=log_file, err_file=err_file)

      attr_row <- inst$attrition
      attr_row$post_outcome_extraction_count   <- harm_res$attrition$post_outcome_extraction_count
      attr_row$post_harmonisation_row_count    <- harm_res$attrition$post_harmonisation_row_count
      attr_row$post_harmonisation_usable_count <- harm_res$attrition$post_harmonisation_usable_count
      attr_row$palindromic_loss_count          <- harm_res$attrition$palindromic_loss_count

      if (is.na(attr_row$non_estimable_reason) && !is.na(harm_res$attrition$non_estimable_reason))
        attr_row$non_estimable_reason <- harm_res$attrition$non_estimable_reason

      attr_row$outcome_name <- o$name
      attr_row$outcome_id   <- o$id
      attr_row$role         <- o$role
      attr_row$phenotype    <- o$phenotype %||% NA_character_
      attr_row$proxy_layer  <- "eQTL"
      attr_row$estimable    <- !is.null(harm_res$harm) &&
        nrow(harm_res$harm) > 0 && sum(harm_res$harm$mr_keep, na.rm=TRUE) > 0

      if (is.na(attr_row$non_estimable_reason) && !isTRUE(attr_row$estimable))
        attr_row$non_estimable_reason <- "insufficient_variant_count_for_estimator"

      all_attr[[length(all_attr)+1]] <- attr_row

      if (isTRUE(attr_row$estimable)) {
        mr_suite <- run_mr_suite(harm=harm_res$harm, outcome_obj=o,
          run_cfg=run_cfg, protocol_cfg=protocol_cfg, log_file=log_file, err_file=err_file)
        if (!is.null(mr_suite$mr) && nrow(mr_suite$mr) > 0) {
          p <- mr_suite$mr; p$proxy_layer <- "eQTL"; all_mr[[length(all_mr)+1]] <- p }
        if (!is.null(mr_suite$heterogeneity) && nrow(mr_suite$heterogeneity) > 0) {
          p <- mr_suite$heterogeneity; p$proxy_layer <- "eQTL"; all_het[[length(all_het)+1]] <- p }
        if (!is.null(mr_suite$pleiotropy) && nrow(mr_suite$pleiotropy) > 0) {
          p <- mr_suite$pleiotropy; p$proxy_layer <- "eQTL"; all_ple[[length(all_ple)+1]] <- p }
        if (!is.null(mr_suite$loo) && nrow(mr_suite$loo) > 0) {
          p <- mr_suite$loo; p$proxy_layer <- "eQTL"; all_loo[[length(all_loo)+1]] <- p }
      }

      log_message(log_file,
        "GENE done gene={g} outcome={o$name} estimable={attr_row$estimable} reason={attr_row$non_estimable_reason %||% 'NA'}")
    }
  }

  list(
    attrition     = dplyr::bind_rows(all_attr),
    mr            = if (length(all_mr)  > 0) dplyr::bind_rows(all_mr)  else data.frame(),
    heterogeneity = if (length(all_het) > 0) dplyr::bind_rows(all_het) else data.frame(),
    pleiotropy    = if (length(all_ple) > 0) dplyr::bind_rows(all_ple) else data.frame(),
    loo           = if (length(all_loo) > 0) dplyr::bind_rows(all_loo) else data.frame()
  )
}

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
write_run_metadata <- function(dirs, protocol_cfg, run_cfg, targets_cfg) {
  write_yaml_snapshot(protocol_cfg, fs::path(dirs$run_dir, "protocol_snapshot.yaml"))
  write_yaml_snapshot(run_cfg,      fs::path(dirs$run_dir, "run_profile_snapshot.yaml"))
  write_yaml_snapshot(targets_cfg,  fs::path(dirs$run_dir, "targets_snapshot.yaml"))
  save_session_info(fs::path(dirs$run_dir, "sessionInfo.txt"))
}

derive_egger_i2gx_table <- function(ple_df) {
  if (is.null(ple_df) || !is.data.frame(ple_df) || nrow(ple_df) == 0) return(data.frame())
  keep <- intersect(c("gene","outcome_name","outcome_id","role","ancestry",
    "egger_intercept","se","pval","i2_gx","egger_reliability",
    "egger_minimum_nsnp","nsnp_usable"), names(ple_df))
  if (length(keep) == 0) return(data.frame())
  ple_df[, keep, drop=FALSE]
}

derive_susie_followup_summary <- function(coloc_df) {
  if (is.null(coloc_df) || !is.data.frame(coloc_df) || nrow(coloc_df) == 0) return(data.frame())
  keep <- intersect(c("gene","outcome_name","method","follow_up_method",
    "susie_triggered","susie_trigger_reason","susie_attempted","susie_status",
    "susie_signal_pairs","susie_PP.H3_max","susie_PP.H4_max",
    "ld_population","ld_reference_panel","ld_source","ld_matrix_status","ld_matrix_note"),
    names(coloc_df))
  if (length(keep) == 0) return(data.frame())
  coloc_df[, keep, drop=FALSE]
}

derive_susie_trigger_log <- function(coloc_df) {
  if (is.null(coloc_df) || !is.data.frame(coloc_df) || nrow(coloc_df) == 0) return(data.frame())
  keep <- intersect(c("gene","outcome_name","susie_triggered","susie_trigger_reason",
    "susie_attempted","susie_status","ld_source","ld_matrix_status","ld_matrix_note"),
    names(coloc_df))
  if (length(keep) == 0) return(data.frame())
  coloc_df[, keep, drop=FALSE]
}

write_core_outputs <- function(dirs, attrition_df, non_estimable_df, mr_df,
                               preferred_df, het_df, ple_df, loo_df,
                               sero_df, transport_df, class_df, class_support_df,
                               priority_coloc_df, overlap_dt=NULL, overlap_summary=NULL) {
  r <- dirs$results_dir
  write_csv_safe(attrition_df,                         fs::path(r, "instrument_attrition_table.csv"))
  write_csv_safe(non_estimable_df,                     fs::path(r, "non_estimable_table.csv"))
  write_csv_safe(mr_df,                                fs::path(r, "main_mr_results.csv"))
  write_csv_safe(preferred_df,                         fs::path(r, "preferred_mr_results.csv"))
  write_csv_safe(het_df,                               fs::path(r, "heterogeneity_results.csv"))
  write_csv_safe(ple_df,                               fs::path(r, "pleiotropy_results.csv"))
  write_csv_safe(derive_egger_i2gx_table(ple_df),      fs::path(r, "egger_i2gx_results.csv"))
  write_csv_safe(loo_df,                               fs::path(r, "leaveoneout_results.csv"))
  write_csv_safe(sero_df,                              fs::path(r, "serostatus_heterogeneity_table.csv"))
  write_csv_safe(transport_df,                         fs::path(r, "ancestry_transportability_table.csv"))
  write_csv_safe(class_df,                             fs::path(r, "benchmark_classification_table.csv"))
  write_csv_safe(class_support_df,                     fs::path(r, "benchmark_classification_with_support_layers.csv"))
  write_csv_safe(priority_coloc_df,                    fs::path(r, "priority_target_colocalisation_summary.csv"))
  write_csv_safe(derive_susie_followup_summary(priority_coloc_df), fs::path(r, "susie_followup_summary.csv"))
  write_csv_safe(derive_susie_trigger_log(priority_coloc_df),      fs::path(r, "susie_trigger_log.csv"))
  if (!is.null(overlap_dt))      write_csv_safe(overlap_dt,      fs::path(r, "sample_overlap_assessment.csv"))
  if (!is.null(overlap_summary)) write_csv_safe(overlap_summary, fs::path(r, "sample_overlap_summary.csv"))
}

write_pqtl_outputs <- function(dirs, pqtl_source_summary, pqtl_res,
                               pqtl_support_df, cross_proxy_df, master_table) {
  safe_df <- function(x) if (is.null(x) || !is.data.frame(x)) data.frame() else x
  r <- dirs$results_dir

  if (is.null(pqtl_res) || !is.list(pqtl_res))
    pqtl_res <- list(attrition=data.frame(), mr=data.frame(), preferred=data.frame(),
                     heterogeneity=data.frame(), pleiotropy=data.frame(), loo=data.frame())

  pqtl_non_est <- if (nrow(safe_df(pqtl_res$attrition)) > 0 &&
                      "estimable" %in% names(pqtl_res$attrition))
    dplyr::filter(pqtl_res$attrition, !estimable) else data.frame()

  write_csv_safe(safe_df(pqtl_source_summary),              fs::path(r, "pqtl_source_summary.csv"))
  write_csv_safe(safe_df(pqtl_res$attrition),               fs::path(r, "pqtl_instrument_attrition_table.csv"))
  write_csv_safe(safe_df(pqtl_res$mr),                      fs::path(r, "pqtl_main_mr_results.csv"))
  write_csv_safe(safe_df(pqtl_res$preferred),               fs::path(r, "pqtl_preferred_mr_results.csv"))
  write_csv_safe(safe_df(pqtl_res$heterogeneity),           fs::path(r, "pqtl_heterogeneity_results.csv"))
  write_csv_safe(safe_df(pqtl_res$pleiotropy),              fs::path(r, "pqtl_pleiotropy_results.csv"))
  write_csv_safe(derive_egger_i2gx_table(safe_df(pqtl_res$pleiotropy)), fs::path(r, "pqtl_egger_i2gx_results.csv"))
  write_csv_safe(safe_df(pqtl_res$loo),                     fs::path(r, "pqtl_leaveoneout_results.csv"))
  write_csv_safe(pqtl_non_est,                              fs::path(r, "pqtl_non_estimable_table.csv"))
  write_csv_safe(safe_df(pqtl_support_df),                  fs::path(r, "pqtl_support_summary.csv"))
  write_csv_safe(safe_df(cross_proxy_df),                   fs::path(r, "cross_proxy_consistency_table.csv"))
  write_csv_safe(safe_df(master_table),                     fs::path(r, "master_benchmark_table.csv"))
}

format_elapsed_seconds <- function(start_time, end_time=Sys.time()) {
  secs <- as.numeric(difftime(end_time, start_time, units="secs"))
  sprintf("%.1f sec (%.2f min)", secs, secs/60)
}

log_stage_timing <- function(log_file, stage_name, start_time, end_time=Sys.time()) {
  log_message(log_file, "STAGE DONE | stage={stage_name} elapsed={format_elapsed_seconds(start_time, end_time)}")
}

# -----------------------------------------------------------------------------
# Main execution function
# -----------------------------------------------------------------------------
main <- function() {
  args          <- commandArgs(trailingOnly=TRUE)
  protocol_path <- if (length(args) >= 1) args[1] else "configs/protocol_v1.yaml"
  run_path      <- if (length(args) >= 2) args[2] else "configs/run_full_study_v1.yaml"

  protocol_cfg <- read_yaml_checked(protocol_path, "Protocol file")
  run_cfg      <- read_yaml_checked(run_path, "Run profile file")
  validate_protocol(protocol_cfg)

  targets_cfg <- read_yaml_checked(protocol_cfg$targets$metadata_file, "Targets metadata file")
  validate_targets_metadata(targets_cfg, protocol_cfg)
  validate_run_profile(run_cfg)
  validate_protocol_run_consistency(protocol_cfg, run_cfg)

  dirs           <- make_run_dirs(prefix=run_cfg$run$run_profile_name)
  run_start_time <- Sys.time()
  log_message(dirs$log_file, "START run_profile={run_cfg$run$run_profile_name}")
  write_run_metadata(dirs, protocol_cfg, run_cfg, targets_cfg)

  target_meta <- read_target_metadata(targets_cfg)
  outcomes    <- flatten_outcomes_for_run(protocol_cfg, run_cfg)
  genes       <- tolower(protocol_cfg$targets$gene_list)
  message(sprintf("Loaded %d genes, %d outcomes from config", length(genes), length(outcomes)))

  # --- Sample overlap ---------------------------------------------------------
  overlap_dt <- NULL; overlap_summary <- NULL
  t0 <- Sys.time()
  if (isTRUE(protocol_cfg$sample_overlap$check_required)) {
    overlap_meta         <- read_yaml_checked(protocol_cfg$sample_overlap$metadata_file, "Sample overlap metadata")
    overlap_exposure_ids <- names(overlap_meta$exposures %||% list())
    if (length(overlap_exposure_ids) == 0)
      overlap_exposure_ids <- overlap_meta$default_exposure_id %||% character(0)
    overlap_dt <- assess_sample_overlap_all_exposures(
      outcomes=outcomes, metadata_file=protocol_cfg$sample_overlap$metadata_file,
      exposure_ids=overlap_exposure_ids)
    overlap_summary <- overlap_dt %>%
      dplyr::count(exposure_id, exposure_name, outcome_name, overlap_status, caveat_level, name="n_rows")
    log_message(dirs$log_file, "Sample overlap completed for exposures={paste(overlap_exposure_ids, collapse=';')}")
  }
  log_stage_timing(dirs$log_file, "sample_overlap", t0)

  # --- Read exposure ----------------------------------------------------------
  t0 <- Sys.time()
  exposure_path <- expand_env_path(protocol_cfg$proxy_layers$primary$exposure_file)
  if (!file.exists(exposure_path))
    stop("Exposure file not found: ", exposure_path,
         "\nCheck INPUT_DATA_PATH env var or protocol_v1.yaml exposure_file setting.")
  exp_raw <- read_exposure_data_standardized(exposure_path)
  write_csv_safe(make_eqtl_gene_sample_summary(exp_raw), fs::path(dirs$results_dir, "eqtl_gene_sample_summary.csv"))
  log_stage_timing(dirs$log_file, "read_exposure", t0)

  # --- eQTL primary suite -----------------------------------------------------
  t0 <- Sys.time()
  eqtl_res <- run_eqtl_primary_suite(
    exp_raw=exp_raw, outcomes=outcomes, genes=genes,
    protocol_cfg=protocol_cfg, run_cfg=run_cfg, target_meta=target_meta,
    log_file=dirs$log_file, err_file=dirs$err_file)
  log_stage_timing(dirs$log_file, "eqtl_primary_suite", t0)

  # --- Preferred / classification ---------------------------------------------
  t0 <- Sys.time()
  preferred_df     <- build_preferred_results(eqtl_res$mr, run_cfg)
  non_estimable_df <- eqtl_res$attrition %>% dplyr::filter(!estimable)
  sero_df          <- make_serostatus_table(preferred_df)
  transport_df     <- make_transportability_table(preferred_df)
  class_df         <- make_benchmark_classification(
    preferred_df=preferred_df, attrition_df=eqtl_res$attrition,
    run_cfg=run_cfg, target_meta=target_meta)
  log_stage_timing(dirs$log_file, "preferred_and_classification", t0)

  # Initialise support-layer placeholders
  priority_coloc_df   <- data.frame()
  pqtl_source_summary <- data.frame()
  pqtl_res <- list(attrition=data.frame(), mr=data.frame(), preferred=data.frame(),
                   heterogeneity=data.frame(), pleiotropy=data.frame(), loo=data.frame())
  pqtl_support_df <- data.frame()
  cross_proxy_df  <- build_cross_proxy_consistency_table(preferred_df, data.frame())

  # --- pQTL validation --------------------------------------------------------
  t0 <- Sys.time()
  if (isTRUE(run_cfg$execution_scope$use_secondary_pqtl_layer)) {
    pqtl_manifest <- read_pqtl_manifest(
      path=run_cfg$pqtl_run$source_manifest, strict=isTRUE(run_cfg$pqtl_run$strict_source_manifest))
    validate_pqtl_manifest(pqtl_manifest, target_genes=genes)
    pqtl_sources        <- load_all_pqtl_sources(pqtl_manifest, target_genes=genes)
    pqtl_source_summary <- make_pqtl_source_summary(pqtl_sources)
    pqtl_res <- run_pqtl_validation_suite(
      pqtl_sources=pqtl_sources, outcomes=outcomes, genes=genes,
      protocol_cfg=protocol_cfg, run_cfg=run_cfg, target_meta=target_meta,
      log_file=dirs$log_file, err_file=dirs$err_file)
    pqtl_support_df <- summarise_pqtl_support(
      preferred_pqtl_df=pqtl_res$preferred, target_meta=target_meta)
    cross_proxy_df <- build_cross_proxy_consistency_table(
      eqtl_preferred_df=preferred_df, pqtl_preferred_df=pqtl_res$preferred)
  }
  log_stage_timing(dirs$log_file, "pqtl_stage", t0)

  # --- Support layers + master table (backbone) -------------------------------
  t0 <- Sys.time()
  class_support_df <- augment_classification_support_layers(
    class_df=class_df, target_meta=target_meta, protocol_cfg=protocol_cfg,
    coloc_df=priority_coloc_df, pqtl_support_df=pqtl_support_df, cross_proxy_df=cross_proxy_df)
  master_table <- build_master_benchmark_table(
    target_meta=target_meta, attrition_df=eqtl_res$attrition,
    eqtl_preferred_df=preferred_df, sero_df=sero_df, transport_df=transport_df,
    class_support_df=class_support_df, coloc_df=priority_coloc_df,
    pqtl_support_df=pqtl_support_df, cross_proxy_df=cross_proxy_df,
    pqtl_attrition_df=pqtl_res$attrition, pqtl_preferred_df=pqtl_res$preferred,
    overlap_dt=overlap_dt)

  # Write backbone outputs BEFORE optional coloc stage
  write_core_outputs(dirs=dirs, attrition_df=eqtl_res$attrition,
    non_estimable_df=non_estimable_df, mr_df=eqtl_res$mr,
    preferred_df=preferred_df, het_df=eqtl_res$heterogeneity,
    ple_df=eqtl_res$pleiotropy, loo_df=eqtl_res$loo,
    sero_df=sero_df, transport_df=transport_df, class_df=class_df,
    class_support_df=class_support_df, priority_coloc_df=priority_coloc_df,
    overlap_dt=overlap_dt, overlap_summary=overlap_summary)
  log_stage_timing(dirs$log_file, "write_backbone_outputs", t0)

  # --- Bias bounds ------------------------------------------------------------
  t0 <- Sys.time()
  eqtl_bias <- if (!is.null(preferred_df) && nrow(preferred_df) > 0)
    dplyr::mutate(preferred_df, exposure_id=protocol_cfg$sample_overlap$default_exposure_id %||% "EQTLGEN_BLOOD")
  else data.frame()
  pqtl_bias <- if (!is.null(pqtl_res$preferred) && nrow(pqtl_res$preferred) > 0)
    dplyr::mutate(pqtl_res$preferred, exposure_id=dplyr::coalesce(
      as.character(.data$pqtl_source_id), as.character(.data$pqtl_source)))
  else data.frame()
  bias_input <- dplyr::bind_rows(eqtl_bias, pqtl_bias)
  if (nrow(bias_input) > 0) {
    write_csv_safe(compute_overlap_bias_bound(bias_input, overlap_fraction=1.0),
      fs::path(dirs$results_dir, "sample_overlap_bias_bound_worst_case.csv"))
    write_csv_safe(compute_overlap_bias_bound(bias_input, overlap_fraction=0.5),
      fs::path(dirs$results_dir, "sample_overlap_bias_bound_50pct.csv"))
  }
  log_stage_timing(dirs$log_file, "bias_bounds", t0)

  # --- pQTL outputs -----------------------------------------------------------
  t0 <- Sys.time()
  write_pqtl_outputs(dirs=dirs, pqtl_source_summary=pqtl_source_summary,
    pqtl_res=pqtl_res, pqtl_support_df=pqtl_support_df,
    cross_proxy_df=cross_proxy_df, master_table=master_table)
  log_stage_timing(dirs$log_file, "write_pqtl_outputs", t0)

  # --- Colocalisation ---------------------------------------------------------
  t0 <- Sys.time()
  if (isTRUE(run_cfg$execution_scope$run_colocalisation)) {
    priority_coloc_df <- run_priority_target_coloc_summary(
      exp_raw=exp_raw, protocol_cfg=protocol_cfg, run_cfg=run_cfg,
      outcomes=outcomes, preferred_df=preferred_df,
      err_file=dirs$err_file, log_file=dirs$log_file)

    class_support_df <- augment_classification_support_layers(
      class_df=class_df, target_meta=target_meta, protocol_cfg=protocol_cfg,
      coloc_df=priority_coloc_df, pqtl_support_df=pqtl_support_df, cross_proxy_df=cross_proxy_df)
    master_table <- build_master_benchmark_table(
      target_meta=target_meta, attrition_df=eqtl_res$attrition,
      eqtl_preferred_df=preferred_df, sero_df=sero_df, transport_df=transport_df,
      class_support_df=class_support_df, coloc_df=priority_coloc_df,
      pqtl_support_df=pqtl_support_df, cross_proxy_df=cross_proxy_df,
      pqtl_attrition_df=pqtl_res$attrition, pqtl_preferred_df=pqtl_res$preferred,
      overlap_dt=overlap_dt)
    master_table <- tri_write_outputs(dirs=dirs, master_table=master_table,
      pqtl_preferred_df=pqtl_res$preferred, overlap_dt=overlap_dt, write_master_augmented=FALSE)

    r <- dirs$results_dir
    write_csv_safe(priority_coloc_df, fs::path(r, "priority_target_colocalisation_summary.csv"))
    write_csv_safe(build_replication_coloc_summary(
        coloc_df=priority_coloc_df, replication_outcome_name="RA_EUR_REP1_FinnGen"),
      fs::path(r, "replication_coloc_summary.csv"))
    write_csv_safe(tri_coloc_cross_outcome_concordance(coloc_df=priority_coloc_df),
      fs::path(r, "coloc_cross_outcome_concordance.csv"))
    write_csv_safe(derive_susie_followup_summary(priority_coloc_df), fs::path(r, "susie_followup_summary.csv"))
    write_csv_safe(derive_susie_trigger_log(priority_coloc_df),      fs::path(r, "susie_trigger_log.csv"))
    write_csv_safe(class_support_df, fs::path(r, "benchmark_classification_with_support_layers.csv"))
    write_csv_safe(master_table,     fs::path(r, "master_benchmark_table.csv"))
  }
  log_stage_timing(dirs$log_file, "colocalisation", t0)

  # --- Failure mode diagnostic ------------------------------------------------
  t0 <- Sys.time()
  if (isTRUE(run_cfg$execution_scope$run_failure_mode_diagnostic)) {
    run_failure_mode_diagnostic(
      results_dir   = dirs$results_dir,
      config_path   = protocol_cfg$failure_mode_diagnostic$config_file %||% "configs/failure_mode_diagnostic.yaml",
      output_dir    = dirs$results_dir,
      log_file      = dirs$log_file,
      err_file      = dirs$err_file,
      run_api_fetch = isTRUE(run_cfg$execution_scope$run_fail_mode_api_fetch))
  }
  log_stage_timing(dirs$log_file, "failure_mode_diagnostic", t0)

  # --- Figures ----------------------------------------------------------------
  write_additional_triangulation_figures(dirs=dirs, master_table=master_table,
    pqtl_source_summary=pqtl_source_summary, priority_coloc_df=priority_coloc_df,
    overlap_dt=overlap_dt, log_file=dirs$log_file)

  t0 <- Sys.time()
  write_full_study_figures(dirs=dirs, protocol_cfg=protocol_cfg,
    attrition_df=eqtl_res$attrition, preferred_df=preferred_df,
    class_support_df=class_support_df, cross_proxy_df=cross_proxy_df, log_file=dirs$log_file)
  log_stage_timing(dirs$log_file, "figures", t0)

  # --- Run manifest -----------------------------------------------------------
  write_csv_safe(
    data.frame(key=c("run_profile_name","protocol_version","purpose","run_dir"),
               value=c(run_cfg$run$run_profile_name, protocol_cfg$project$protocol_version,
                       run_cfg$run$purpose, as.character(dirs$run_dir)),
               stringsAsFactors=FALSE),
    fs::path(dirs$run_dir, "run_manifest.csv"))

  log_stage_timing(dirs$log_file, "entire_run", run_start_time)
  log_message(dirs$log_file, "DONE run_dir={dirs$run_dir}")
  cat("Run completed:", as.character(dirs$run_dir), "\n")
}

# -----------------------------------------------------------------------------
# Entry point
# -----------------------------------------------------------------------------
main()
