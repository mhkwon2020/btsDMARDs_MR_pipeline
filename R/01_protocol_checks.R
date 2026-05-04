# R/01_protocol_checks.R
# Protocol validation utilities
# Version: 1.0
# Author: Mihye Kwon
suppressPackageStartupMessages({
  library(yaml)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

fail <- function(fmt, ...) {
  stop(sprintf(fmt, ...), call. = FALSE)
}

check_nonempty <- function(x, label) {
  if (is.null(x)) {
    fail("%s is missing", label)
  }
  if (is.character(x) && length(x) == 1 && !nzchar(trimws(x))) {
    fail("%s must not be empty", label)
  }
  if (length(x) == 0) {
    fail("%s must not be empty", label)
  }
  invisible(TRUE)
}

check_true_false <- function(x, label) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    fail("%s must be a single TRUE/FALSE value", label)
  }
  invisible(TRUE)
}

check_scalar_numeric <- function(x, label, allow_na = FALSE) {
  if (allow_na && (is.null(x) || length(x) == 0 || is.na(x))) return(invisible(TRUE))
  if (!is.numeric(x) || length(x) != 1 || is.na(x)) {
    fail("%s must be a single numeric value", label)
  }
  invisible(TRUE)
}

check_file_exists <- function(path, label) {
  check_nonempty(path, label)
  if (!file.exists(path)) {
    fail("%s does not exist: %s", label, path)
  }
  invisible(TRUE)
}

read_yaml_checked <- function(path, label = "YAML file") {
  check_file_exists(path, label)
  out <- tryCatch(
    yaml::read_yaml(path),
    error = function(e) {
      fail("%s could not be parsed as YAML: %s", label, e$message)
    }
  )
  if (is.null(out)) {
    fail("%s could not be parsed as YAML: %s", label, path)
  }
  out
}

validate_protocol <- function(cfg) {
  required_sections <- c(
    "project", "design", "targets", "proxy_layers",
    "instrument_selection", "outcomes", "harmonisation",
    "mr_estimation", "assumption_framework", "sensitivity_analysis",
    "colocalisation", "non_estimable_definition",
    "benchmark_classification", "reporting", "environment"
  )

  for (sec in required_sections) {
    check_nonempty(cfg[[sec]], paste0("protocol.", sec))
  }

  check_nonempty(cfg$project$name, "protocol.project.name")
  check_nonempty(cfg$project$protocol_version, "protocol.project.protocol_version")
  check_nonempty(cfg$project$study_identity, "protocol.project.study_identity")

  check_true_false(cfg$targets$prespecified, "protocol.targets.prespecified")
  check_true_false(cfg$targets$allow_dynamic_target_addition,
                   "protocol.targets.allow_dynamic_target_addition")
  check_nonempty(cfg$targets$metadata_file, "protocol.targets.metadata_file")
  check_scalar_numeric(cfg$targets$expected_target_count, "protocol.targets.expected_target_count")
  check_nonempty(cfg$targets$gene_list, "protocol.targets.gene_list")
  check_file_exists(cfg$targets$metadata_file, "Target metadata file")

  check_true_false(cfg$proxy_layers$primary$enabled, "protocol.proxy_layers.primary.enabled")
  check_nonempty(cfg$proxy_layers$primary$type, "protocol.proxy_layers.primary.type")
  check_nonempty(cfg$proxy_layers$primary$role, "protocol.proxy_layers.primary.role")
  check_nonempty(cfg$proxy_layers$primary$exposure_file, "protocol.proxy_layers.primary.exposure_file")

  if (!is.null(cfg$proxy_layers$secondary)) {
    check_true_false(cfg$proxy_layers$secondary$enabled,
                     "protocol.proxy_layers.secondary.enabled")
    check_nonempty(cfg$proxy_layers$secondary$type,
                   "protocol.proxy_layers.secondary.type")
    check_nonempty(cfg$proxy_layers$secondary$role,
                   "protocol.proxy_layers.secondary.role")
  }

  check_nonempty(cfg$instrument_selection$p_threshold, "protocol.instrument_selection.p_threshold")
  check_scalar_numeric(cfg$instrument_selection$p_threshold$primary,
                       "protocol.instrument_selection.p_threshold.primary")
  check_nonempty(cfg$instrument_selection$clumping, "protocol.instrument_selection.clumping")
  check_scalar_numeric(cfg$instrument_selection$clumping$r2,
                       "protocol.instrument_selection.clumping.r2")
  check_scalar_numeric(cfg$instrument_selection$clumping$kb,
                       "protocol.instrument_selection.clumping.kb")

  check_nonempty(cfg$outcomes$primary, "protocol.outcomes.primary")
  check_nonempty(cfg$outcomes$replication, "protocol.outcomes.replication")
  check_nonempty(cfg$outcomes$secondary, "protocol.outcomes.secondary")

  check_true_false(cfg$sample_overlap$check_required, "protocol.sample_overlap.check_required")
  check_nonempty(cfg$sample_overlap$metadata_file, "protocol.sample_overlap.metadata_file")
  check_file_exists(cfg$sample_overlap$metadata_file, "Sample overlap metadata file")

  check_scalar_numeric(cfg$harmonisation$action, "protocol.harmonisation.action")
  check_true_false(cfg$mr_estimation$rule_based_method_selection,
                 "protocol.mr_estimation.rule_based_method_selection")
  check_nonempty(cfg$mr_estimation$primary_estimator,
                 "protocol.mr_estimation.primary_estimator")
  check_nonempty(cfg$mr_estimation$method_rules,
                 "protocol.mr_estimation.method_rules")
  check_nonempty(cfg$mr_estimation$sensitivity_estimators,
                 "protocol.mr_estimation.sensitivity_estimators")

  check_true_false(cfg$sensitivity_analysis$threshold_sensitivity,
                   "protocol.sensitivity_analysis.threshold_sensitivity")
  check_true_false(cfg$sensitivity_analysis$estimator_sensitivity,
                   "protocol.sensitivity_analysis.estimator_sensitivity")
  check_true_false(cfg$sensitivity_analysis$cross_dataset_replication,
                   "protocol.sensitivity_analysis.cross_dataset_replication")
  check_true_false(cfg$sensitivity_analysis$serostatus_heterogeneity,
                   "protocol.sensitivity_analysis.serostatus_heterogeneity")
  check_true_false(cfg$sensitivity_analysis$ancestry_transportability,
                   "protocol.sensitivity_analysis.ancestry_transportability")
  check_true_false(cfg$sensitivity_analysis$pQTL_validation,
                   "protocol.sensitivity_analysis.pQTL_validation")
  check_true_false(cfg$sensitivity_analysis$colocalisation_support,
                   "protocol.sensitivity_analysis.colocalisation_support")

  check_true_false(cfg$colocalisation$enabled, "protocol.colocalisation.enabled")
  check_nonempty(cfg$colocalisation$priority_targets, "protocol.colocalisation.priority_targets")
  if (!is.null(cfg$colocalisation$susie_follow_up)) {
    check_true_false(cfg$colocalisation$susie_follow_up$enabled,
                     "protocol.colocalisation.susie_follow_up.enabled")
    check_scalar_numeric(cfg$colocalisation$susie_follow_up$minimum_overlap_snps,
                         "protocol.colocalisation.susie_follow_up.minimum_overlap_snps")
    check_scalar_numeric(cfg$colocalisation$susie_follow_up$minimum_ld_snps,
                         "protocol.colocalisation.susie_follow_up.minimum_ld_snps")
  }
  if (!is.null(cfg$colocalisation$ld_matrix)) {
    check_nonempty(cfg$colocalisation$ld_matrix$source_priority,
                   "protocol.colocalisation.ld_matrix.source_priority")
    local_cfg <- cfg$colocalisation$ld_matrix$local_plink
    if (is.null(local_cfg)) {
      stop("protocol.colocalisation.ld_matrix.local_plink is required for local LD")
    }
    check_true_false(local_cfg$enabled,
                     "protocol.colocalisation.ld_matrix.local_plink.enabled")
    check_nonempty(local_cfg$plink_bin,
                   "protocol.colocalisation.ld_matrix.local_plink.plink_bin")
    check_nonempty(local_cfg$reference_bfiles$EUR,
                   "protocol.colocalisation.ld_matrix.local_plink.reference_bfiles.EUR")
    check_nonempty(local_cfg$reference_bfiles$EAS,
                   "protocol.colocalisation.ld_matrix.local_plink.reference_bfiles.EAS")
    src_priority <- unlist(cfg$colocalisation$ld_matrix$source_priority)
    if (!identical(as.character(src_priority), "local_plink")) {
      stop("protocol.colocalisation.ld_matrix.source_priority must be local_plink only")
    }
  }

  check_true_false(cfg$non_estimable_definition$treat_as_negative_evidence,
                   "protocol.non_estimable_definition.treat_as_negative_evidence")
  check_nonempty(cfg$benchmark_classification$classes,
                 "protocol.benchmark_classification.classes")
  check_nonempty(cfg$benchmark_classification$dimensions,
                 "protocol.benchmark_classification.dimensions")

  invisible(TRUE)
}

validate_targets_metadata <- function(targets_cfg, protocol_cfg = NULL) {
  check_nonempty(targets_cfg$version, "targets.version")
  check_nonempty(targets_cfg$targets, "targets.targets")

  target_genes <- vapply(targets_cfg$targets, function(x) tolower(x$gene %||% NA_character_), character(1))
  if (anyNA(target_genes) || any(!nzchar(target_genes))) {
    fail("Every target entry in targets.yaml must define a non-empty gene field")
  }
  if (anyDuplicated(target_genes)) {
    fail("targets.yaml contains duplicated gene entries: %s",
         paste(unique(target_genes[duplicated(target_genes)]), collapse = ", "))
  }

  if (!is.null(protocol_cfg)) {
    protocol_genes <- tolower(protocol_cfg$targets$gene_list)
    if (!setequal(protocol_genes, target_genes)) {
      missing_in_targets <- setdiff(protocol_genes, target_genes)
      extra_in_targets <- setdiff(target_genes, protocol_genes)
      fail(
        paste0(
          "Protocol target list and targets metadata do not match. ",
          "Missing in targets metadata: %s. Extra in targets metadata: %s."
        ),
        paste(missing_in_targets, collapse = ", ") %||% "none",
        paste(extra_in_targets, collapse = ", ") %||% "none"
      )
    }
    if (length(protocol_genes) != protocol_cfg$targets$expected_target_count) {
      fail("protocol.targets.expected_target_count (%s) does not match protocol.targets.gene_list length (%s)",
           protocol_cfg$targets$expected_target_count,
           length(protocol_genes))
    }
  }

  invisible(TRUE)
}

get_required_outputs <- function(run_cfg) {
  run_cfg$required_outputs_full_study %||%
    run_cfg$required_outputs_today %||%
    run_cfg$required_outputs %||%
    character(0)
}

validate_run_profile <- function(run_cfg) {
  required_sections <- c(
    "run", "execution_scope", "targets", "thresholds",
    "outcomes_to_run", "mr_methods", "sensitivity_to_run",
    "colocalisation_run", "reporting"
  )
  for (sec in required_sections) {
    check_nonempty(run_cfg[[sec]], paste0("run.", sec))
  }

  check_nonempty(run_cfg$run$run_profile_name, "run.run.run_profile_name")
  check_nonempty(run_cfg$run$protocol_version, "run.run.protocol_version")
  check_nonempty(run_cfg$run$purpose, "run.run.purpose")

  check_true_false(run_cfg$execution_scope$use_primary_eqtl_layer,
                   "run.execution_scope.use_primary_eqtl_layer")
  check_true_false(run_cfg$execution_scope$use_secondary_pqtl_layer,
                   "run.execution_scope.use_secondary_pqtl_layer")
  check_true_false(run_cfg$execution_scope$run_colocalisation,
                   "run.execution_scope.run_colocalisation")
  check_nonempty(run_cfg$execution_scope$colocalisation_mode,
                 "run.execution_scope.colocalisation_mode")

  valid_coloc_modes <- c("priority_targets_primary_outcome",
                         "priority_targets_multi_outcome")
  if (!run_cfg$execution_scope$colocalisation_mode %in% valid_coloc_modes) {
    fail("run.execution_scope.colocalisation_mode must be one of: %s",
         paste(valid_coloc_modes, collapse = ", "))
  }
  if (identical(run_cfg$execution_scope$colocalisation_mode,
                "priority_targets_multi_outcome")) {
    check_nonempty(run_cfg$colocalisation_run$coloc_outcomes,
                   "run.colocalisation_run.coloc_outcomes")
    check_nonempty(run_cfg$colocalisation_run$coloc_outcomes$core,
                   "run.colocalisation_run.coloc_outcomes.core")
  }

  check_true_false(run_cfg$targets$use_all_prespecified_targets,
                   "run.targets.use_all_prespecified_targets")

  check_true_false(run_cfg$thresholds$use_primary_p_threshold_only,
                   "run.thresholds.use_primary_p_threshold_only")
  check_true_false(run_cfg$thresholds$run_threshold_sensitivity_today,
                   "run.thresholds.run_threshold_sensitivity_today")

  check_true_false(run_cfg$outcomes_to_run$primary, "run.outcomes_to_run.primary")
  check_true_false(run_cfg$outcomes_to_run$replication, "run.outcomes_to_run.replication")
  check_nonempty(run_cfg$outcomes_to_run$secondary, "run.outcomes_to_run.secondary")
  check_true_false(run_cfg$outcomes_to_run$secondary$seropos,
                   "run.outcomes_to_run.secondary.seropos")
  check_true_false(run_cfg$outcomes_to_run$secondary$seroneg,
                   "run.outcomes_to_run.secondary.seroneg")
  check_true_false(run_cfg$outcomes_to_run$secondary$eas,
                   "run.outcomes_to_run.secondary.eas")

  check_true_false(run_cfg$mr_methods$apply_rule_based_method_selection_from_protocol,
                   "run.mr_methods.apply_rule_based_method_selection_from_protocol")

  check_true_false(run_cfg$sensitivity_to_run$heterogeneity,
                   "run.sensitivity_to_run.heterogeneity")
  check_true_false(run_cfg$sensitivity_to_run$pleiotropy,
                   "run.sensitivity_to_run.pleiotropy")
  check_true_false(run_cfg$sensitivity_to_run$leave_one_out,
                   "run.sensitivity_to_run.leave_one_out")
  check_true_false(run_cfg$sensitivity_to_run$estimator_sensitivity,
                   "run.sensitivity_to_run.estimator_sensitivity")
  check_true_false(run_cfg$sensitivity_to_run$threshold_sensitivity,
                   "run.sensitivity_to_run.threshold_sensitivity")
  if (!is.null(run_cfg$sensitivity_to_run$report_i2_gx_with_egger)) {
    check_true_false(run_cfg$sensitivity_to_run$report_i2_gx_with_egger,
                     "run.sensitivity_to_run.report_i2_gx_with_egger")
  }

  check_true_false(run_cfg$colocalisation_run$enabled, "run.colocalisation_run.enabled")
  check_nonempty(run_cfg$colocalisation_run$targets, "run.colocalisation_run.targets")
  check_nonempty(run_cfg$colocalisation_run$method_default,
                 "run.colocalisation_run.method_default")

  if (!is.null(run_cfg$colocalisation_run$allow_extended_method_today)) {
    check_true_false(run_cfg$colocalisation_run$allow_extended_method_today,
                     "run.colocalisation_run.allow_extended_method_today")
  }
  if (!is.null(run_cfg$colocalisation_run$require_ld_for_susie)) {
    check_true_false(run_cfg$colocalisation_run$require_ld_for_susie,
                     "run.colocalisation_run.require_ld_for_susie")
  }
  if (!is.null(run_cfg$colocalisation_run$ld_backend_for_susie)) {
    if (!identical(run_cfg$colocalisation_run$ld_backend_for_susie, "local_plink_only")) {
      stop("run.colocalisation_run.ld_backend_for_susie must be local_plink_only")
    }
  }

  required_outputs <- get_required_outputs(run_cfg)
  if (length(required_outputs) == 0) {
    fail("Run profile must define at least one required output")
  }

  if (!is.null(run_cfg$abstract_priority_rules)) {
    check_true_false(run_cfg$abstract_priority_rules$require_il6r_non_estimable_explanation,
                     "run.abstract_priority_rules.require_il6r_non_estimable_explanation")
    check_true_false(run_cfg$abstract_priority_rules$require_benchmark_framing,
                     "run.abstract_priority_rules.require_benchmark_framing")
    check_true_false(run_cfg$abstract_priority_rules$require_positive_control_interpretation,
                     "run.abstract_priority_rules.require_positive_control_interpretation")
  }

  if (!is.null(run_cfg$pqtl_run)) {
    check_true_false(run_cfg$pqtl_run$enabled, "run.pqtl_run.enabled")
    if (isTRUE(run_cfg$pqtl_run$enabled)) {
      check_nonempty(run_cfg$pqtl_run$source_manifest, "run.pqtl_run.source_manifest")
      if (!is.null(run_cfg$pqtl_run$strict_source_manifest)) {
        check_true_false(run_cfg$pqtl_run$strict_source_manifest,
                         "run.pqtl_run.strict_source_manifest")
      }
      if (!is.null(run_cfg$pqtl_run$run_replication_outcomes)) {
        check_true_false(run_cfg$pqtl_run$run_replication_outcomes,
                         "run.pqtl_run.run_replication_outcomes")
      }
    }
  }

  invisible(TRUE)
}

validate_protocol_run_consistency <- function(protocol_cfg, run_cfg) {
  if (!identical(protocol_cfg$project$protocol_version, run_cfg$run$protocol_version)) {
    fail("Protocol version mismatch: protocol=%s, run=%s",
         protocol_cfg$project$protocol_version, run_cfg$run$protocol_version)
  }

  if (!isTRUE(run_cfg$execution_scope$use_primary_eqtl_layer)) {
    fail("run.execution_scope.use_primary_eqtl_layer must be TRUE")
  }

  if (!isTRUE(protocol_cfg$proxy_layers$primary$enabled)) {
    fail("protocol.proxy_layers.primary.enabled must be TRUE")
  }

  if (!identical(protocol_cfg$proxy_layers$primary$type, "cis_eQTL")) {
    fail("protocol.proxy_layers.primary.type must be 'cis_eQTL'")
  }

  if (isTRUE(run_cfg$execution_scope$use_secondary_pqtl_layer)) {
    if (!isTRUE(protocol_cfg$proxy_layers$secondary$enabled)) {
      fail("Run requests secondary pQTL layer but protocol secondary layer is not enabled")
    }
    if (!isTRUE(protocol_cfg$sensitivity_analysis$pQTL_validation)) {
      fail("Run requests secondary pQTL layer but protocol.sensitivity_analysis.pQTL_validation is FALSE")
    }
    if (is.null(run_cfg$pqtl_run) || !isTRUE(run_cfg$pqtl_run$enabled)) {
      fail("Run requests secondary pQTL layer but run.pqtl_run.enabled is not TRUE")
    }
  }

  if (isTRUE(run_cfg$execution_scope$run_colocalisation)) {
    if (!isTRUE(protocol_cfg$colocalisation$enabled)) {
      fail("Run requests colocalisation but protocol.colocalisation.enabled is FALSE")
    }
    if (!isTRUE(run_cfg$colocalisation_run$enabled)) {
      fail("Run requests colocalisation but run.colocalisation_run.enabled is FALSE")
    }
  }

  if (!all(run_cfg$colocalisation_run$targets %in% protocol_cfg$colocalisation$priority_targets)) {
    fail("run.colocalisation_run.targets must be a subset of protocol.colocalisation.priority_targets")
  }

  if (isTRUE(run_cfg$thresholds$run_threshold_sensitivity_today) &&
      !isTRUE(protocol_cfg$instrument_selection$p_threshold$sensitivity$enabled)) {
    fail("Run requests threshold sensitivity but protocol threshold sensitivity is disabled")
  }

  if (isTRUE(run_cfg$sensitivity_to_run$threshold_sensitivity) &&
      !isTRUE(protocol_cfg$sensitivity_analysis$threshold_sensitivity)) {
    fail("Run requests threshold sensitivity but protocol.sensitivity_analysis.threshold_sensitivity is FALSE")
  }

  if (!isTRUE(run_cfg$mr_methods$apply_rule_based_method_selection_from_protocol)) {
    fail("run.mr_methods.apply_rule_based_method_selection_from_protocol must be TRUE")
  }

  required_outputs <- get_required_outputs(run_cfg)

  required_core <- c(
    "instrument_attrition_table",
    "non_estimable_table",
    "main_mr_results",
    "preferred_mr_results",
    "benchmark_classification_table"
  )

  missing_core <- setdiff(required_core, required_outputs)
  if (length(missing_core) > 0) {
    fail("Run profile is missing required core outputs: %s",
         paste(missing_core, collapse = ", "))
  }

  if (isTRUE(run_cfg$execution_scope$run_colocalisation) &&
      !"priority_target_colocalisation_summary" %in% required_outputs) {
    fail("Run enables colocalisation but required outputs do not include priority_target_colocalisation_summary")
  }

  if (isTRUE(run_cfg$execution_scope$use_secondary_pqtl_layer)) {
    pqtl_required <- c("pqtl_main_mr_results", "pqtl_preferred_mr_results", "pqtl_support_summary")
    missing_pqtl <- setdiff(pqtl_required, required_outputs)
    if (length(missing_pqtl) > 0) {
      fail("Run enables pQTL but required outputs do not include: %s",
           paste(missing_pqtl, collapse = ", "))
    }
    if ("pqtl_egger_i2gx_results" %in% setdiff(c("pqtl_egger_i2gx_results"), required_outputs)) {
      fail("Run enables pQTL but required outputs do not include pqtl_egger_i2gx_results")
    }
  } else {
    unexpected_pqtl <- grep("^pqtl_", required_outputs, value = TRUE)
    if (length(unexpected_pqtl) > 0) {
      fail("Run disables pQTL but still requires pQTL outputs: %s",
           paste(unexpected_pqtl, collapse = ", "))
    }
  }

  if (isTRUE(run_cfg$sensitivity_to_run$pleiotropy) &&
      isTRUE(run_cfg$sensitivity_to_run$report_i2_gx_with_egger %||% TRUE) &&
      !"egger_i2gx_results" %in% required_outputs) {
    fail("Run enables Egger/I2_GX reporting but required outputs do not include egger_i2gx_results")
  }

  if (isTRUE(run_cfg$execution_scope$run_colocalisation) &&
      isTRUE(run_cfg$colocalisation_run$allow_extended_method_today %||% FALSE)) {
    susie_required <- c("susie_followup_summary", "susie_trigger_log")
    missing_susie <- setdiff(susie_required, required_outputs)
    if (length(missing_susie) > 0) {
      fail("Run enables extended coloc follow-up but required outputs do not include: %s",
           paste(missing_susie, collapse = ", "))
    }
  }

  if (identical(run_cfg$execution_scope$colocalisation_mode,
                "priority_targets_multi_outcome") &&
      length(run_cfg$colocalisation_run$coloc_outcomes$supplementary %||% character(0)) > 0) {
    if (!"replication_coloc_summary" %in% required_outputs) {
      fail("Run enables multi-outcome coloc with supplementary outcomes but required outputs do not include replication_coloc_summary")
    }
  }

  invisible(TRUE)
}

print_summary <- function(protocol_cfg, targets_cfg, run_cfg) {
  cat("Protocol, targets metadata, and run profile validation passed.\n")
  cat(sprintf("Project: %s\n", protocol_cfg$project$name))
  cat(sprintf("Protocol version: %s\n", protocol_cfg$project$protocol_version))
  cat(sprintf("Run profile: %s\n", run_cfg$run$run_profile_name))
  cat(sprintf("Purpose: %s\n", run_cfg$run$purpose))
  cat(sprintf("Primary exposure file: %s\n", protocol_cfg$proxy_layers$primary$exposure_file))
  cat(sprintf("Number of targets in protocol: %d\n", length(protocol_cfg$targets$gene_list)))
  cat(sprintf("Number of targets in metadata: %d\n", length(targets_cfg$targets)))
  cat(sprintf("Primary p-threshold: %s\n",
              as.character(protocol_cfg$instrument_selection$p_threshold$primary)))
  cat(sprintf("Today's run: eQTL=%s, pQTL=%s, colocalisation=%s\n",
              run_cfg$execution_scope$use_primary_eqtl_layer,
              run_cfg$execution_scope$use_secondary_pqtl_layer,
              run_cfg$execution_scope$run_colocalisation))
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  protocol_path <- if (length(args) >= 1) args[1] else "configs/protocol_v1.yaml"
  run_path <- if (length(args) >= 2) args[2] else "configs/run_full_study_v1.yaml"

  protocol_cfg <- read_yaml_checked(protocol_path, "Protocol file")
  run_cfg <- read_yaml_checked(run_path, "Run profile file")

  validate_protocol(protocol_cfg)

  targets_path <- protocol_cfg$targets$metadata_file
  targets_cfg <- read_yaml_checked(targets_path, "Targets metadata file")
  validate_targets_metadata(targets_cfg, protocol_cfg)

  validate_run_profile(run_cfg)
  validate_protocol_run_consistency(protocol_cfg, run_cfg)
  print_summary(protocol_cfg, targets_cfg, run_cfg)
}

if (sys.nframe() == 0) {
  main()
}
