# R/02_input_readers.R
# Input reader utilities
# Version: 1.0
# Author: Mihye Kwon
suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
})

read_protocol_bundle <- function(protocol_path = "configs/protocol_v1.yaml",
                                 run_path = "configs/run_full_study_v1.yaml") {
  protocol_cfg <- yaml::read_yaml(protocol_path)
  run_cfg <- yaml::read_yaml(run_path)
  targets_cfg <- yaml::read_yaml(protocol_cfg$targets$metadata_file)
  list(protocol = protocol_cfg, run = run_cfg, targets = targets_cfg)
}

flatten_outcomes_for_run <- function(protocol_cfg, run_cfg) {
  out <- list()

  add_items <- function(items, role) {
    lapply(items, function(x) {
      x$role <- role
      x
    })
  }

  if (isTRUE(run_cfg$outcomes_to_run$primary)) {
    out <- c(out, add_items(protocol_cfg$outcomes$primary, "primary"))
  }
  if (isTRUE(run_cfg$outcomes_to_run$replication)) {
    out <- c(out, add_items(protocol_cfg$outcomes$replication, "replication"))
  }
  if (isTRUE(run_cfg$outcomes_to_run$secondary$seropos)) {
    out <- c(out, add_items(
      Filter(function(x) identical(x$name, "RA_EUR_SEROPOS"), protocol_cfg$outcomes$secondary),
      "secondary"
    ))
  }
  if (isTRUE(run_cfg$outcomes_to_run$secondary$seroneg)) {
    out <- c(out, add_items(
      Filter(function(x) identical(x$name, "RA_EUR_SERONEG"), protocol_cfg$outcomes$secondary),
      "secondary"
    ))
  }
  if (isTRUE(run_cfg$outcomes_to_run$secondary$eas)) {
    out <- c(out, add_items(
      Filter(function(x) identical(x$name, "RA_EAS_BBJ"), protocol_cfg$outcomes$secondary),
      "secondary"
    ))
  }

  out
}

derive_effect_allele_frequency <- function(df) {
  needed <- c("AssessedAllele", "AlleleA", "AlleleB", "AlleleB_all")
  miss <- setdiff(needed, names(df))
  if (length(miss) > 0) {
    stop("Cannot derive effect allele frequency. Missing columns: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }

  out <- ifelse(
    df$AssessedAllele == df$AlleleB,
    df$AlleleB_all,
    ifelse(df$AssessedAllele == df$AlleleA, 1 - df$AlleleB_all, NA_real_)
  )

  out
}

standardize_exposure_eqtlgen <- function(df) {
  needed <- c(
    "SNP", "Pvalue", "SNPChr", "SNPPos",
    "AssessedAllele", "OtherAllele",
    "GeneSymbol", "Beta", "SE",
    "AlleleA", "AlleleB", "AlleleB_all", "NrSamples"
  )
  miss <- setdiff(needed, names(df))
  if (length(miss) > 0) {
    stop("Exposure file missing required columns: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }

  df2 <- df %>%
    mutate(
      gene = tolower(GeneSymbol),
      beta = as.numeric(Beta),
      se = as.numeric(SE),
      pval = as.numeric(Pvalue),
      chr = as.integer(SNPChr),
      pos = as.numeric(SNPPos),
      effect_allele = AssessedAllele,
      other_allele = OtherAllele,
      n_eqtl_variant = as.numeric(NrSamples)
    )

  df2$eaf <- derive_effect_allele_frequency(df2)
  df2$maf <- pmin(df2$eaf, 1 - df2$eaf)

  gene_n_summary <- df2 %>%
    group_by(gene) %>%
    summarise(
      n_eqtl_gene_median = stats::median(n_eqtl_variant, na.rm = TRUE),
      n_eqtl_gene_max = max(n_eqtl_variant, na.rm = TRUE),
      n_eqtl_gene_min = min(n_eqtl_variant, na.rm = TRUE),
      n_eqtl_gene_distinct = n_distinct(n_eqtl_variant),
      .groups = "drop"
    )

  df2 <- df2 %>%
    left_join(gene_n_summary, by = "gene")

  df2
}

read_exposure_data_standardized <- function(path) {
  raw <- utils::read.table(
    path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  standardize_exposure_eqtlgen(raw)
}

read_target_metadata <- function(targets_cfg) {
  rows <- lapply(targets_cfg$targets, function(x) {
    data.frame(
      gene = tolower(x$gene %||% NA_character_),
      class = x$class %||% NA_character_,
      target_role = x$target_role %||% NA_character_,
      tier = x$tier %||% NA_character_,
      direct_ra_target = x$direct_ra_target %||% NA,
      drug_examples = paste(x$drug_examples %||% character(0), collapse = "; "),
      drug_mechanism = x$drug_mechanism %||% NA_character_,
      expected_protective_direction_if_aligned = x$expected_protective_direction_if_aligned %||% NA_character_,
      alignment_certainty = x$alignment_certainty %||% NA_character_,
      interpretability_note = x$interpretability_note %||% NA_character_,
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(rows)
}

make_eqtl_gene_sample_summary <- function(exp_raw) {
  exp_raw %>%
    distinct(
      gene,
      n_eqtl_gene_median,
      n_eqtl_gene_max,
      n_eqtl_gene_min,
      n_eqtl_gene_distinct
    ) %>%
    arrange(gene)
}
