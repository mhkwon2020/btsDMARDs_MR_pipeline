# R/13_triangulation.R
# Triangulation utilities
# Version: 1.0
# Author: Mihye Kwon
suppressPackageStartupMessages({
  library(dplyr)
})

# -----------------------------------------------------------------------------
# R/13_triangulation.R
# Additive triangulation utilities for the btsDMARDs_MR full-study pipeline.
#
# Design principles:
# - Do NOT override existing functions in R/11_reporting_tables.R
# - Use only additive, uniquely-prefixed helpers (`tri_`)
# - Be robust to minor column-name variation in master_benchmark_table
# - Avoid assuming unavailable columns; degrade gracefully to NA
# -----------------------------------------------------------------------------

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
  }
}

tri_pick_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

tri_get_vec <- function(df, candidates, default = NA) {
  col <- tri_pick_col(df, candidates)
  if (is.na(col)) return(rep(default, nrow(df)))
  df[[col]]
}

tri_first_non_missing <- function(...) {
  vals <- unlist(list(...), use.names = FALSE)
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) return(NA)
  vals[1]
}

tri_collapse_notes <- function(...) {
  vals <- unlist(list(...), use.names = FALSE)
  vals <- trimws(as.character(vals))
  vals <- vals[nzchar(vals) & !is.na(vals)]
  if (length(vals) == 0) return(NA_character_)
  paste(unique(vals), collapse = " | ")
}

tri_safe_bool <- function(x) {
  ifelse(is.na(x), NA, isTRUE(x))
}

tri_choose_coloc_support <- function(master_table) {
  if (!is.data.frame(master_table) || nrow(master_table) == 0) return(character(0))
  col <- tri_pick_col(master_table, c("final_coloc_support", "coloc_support"))
  if (is.na(col)) return(rep(NA_character_, nrow(master_table)))
  as.character(master_table[[col]])
}

tri_choose_gene_col <- function(master_table) {
  tri_pick_col(master_table, c("gene", "target_gene"))
}

tri_choose_eqtl_class_col <- function(master_table) {
  tri_pick_col(master_table, c("class.y", "class", "benchmark_class", "eqtl_class", "mr_class"))
}

tri_classify_summary_class <- function(eqtl_class,
                                       eqtl_primary_support,
                                       eqtl_replication_support,
                                       pqtl_ukb_support,
                                       pqtl_decode_support,
                                       cross_source_robustness,
                                       coloc_support,
                                       discordance_type,
                                       technical_caution) {
  coloc_strong <- !is.na(coloc_support) && coloc_support %in% c("Strong_coloc")
  any_pqtl_support <- isTRUE(pqtl_ukb_support) || isTRUE(pqtl_decode_support)

  if (!is.na(eqtl_class) && identical(eqtl_class, "Non_estimable")) {
    return("Backbone_non_estimable")
  }

  if (!is.na(cross_source_robustness) && identical(cross_source_robustness, "discordant_support_two_source")) {
    return("Cross_source_discordant")
  }

  if (!is.na(discordance_type) && grepl("Type3", discordance_type, fixed = TRUE)) {
    return("Cross_source_discordant")
  }

  if (!is.na(eqtl_class) && identical(eqtl_class, "Strong") && any_pqtl_support && coloc_strong) {
    return("Strong_triangulation")
  }

  if (!is.na(eqtl_class) && identical(eqtl_class, "Strong") && any_pqtl_support) {
    return("Strong_with_partial_orthogonal_support")
  }

  if (!is.na(eqtl_class) && identical(eqtl_class, "Strong") && coloc_strong) {
    return("Strong_eqtl_coloc")
  }

  if (!is.na(eqtl_class) && identical(eqtl_class, "Weak_or_partial") && any_pqtl_support) {
    return("Partial_triangulation")
  }

  if (!is.na(eqtl_class) && eqtl_class %in% c("Strong", "Weak_or_partial")) {
    return("eQTL_led_support")
  }

  if (!is.na(discordance_type) && nzchar(discordance_type)) {
    return("Discordant_evidence_pattern")
  }

  if (isTRUE(technical_caution)) {
    return("Cautious_not_supported")
  }

  "Not_supported_overall"
}

tri_add_summary_class <- function(master_table) {
  if (is.null(master_table) || !is.data.frame(master_table) || nrow(master_table) == 0) {
    return(data.frame())
  }

  eqtl_class <- tri_get_vec(master_table, c("class.y", "class", "benchmark_class", "eqtl_class"), default = NA_character_)
  eqtl_primary_support <- tri_get_vec(master_table, c("primary_nominal_support", "eqtl_primary_support", "primary_support"), default = NA)
  eqtl_replication_support <- tri_get_vec(master_table, c("replication_nominal_support", "eqtl_replication_support", "replication_support"), default = NA)
  pqtl_ukb_support <- tri_get_vec(master_table, c("pqtl_ukbppp_primary_support"), default = NA)
  pqtl_decode_support <- tri_get_vec(master_table, c("pqtl_decode_primary_support"), default = NA)
  cross_source_robustness <- tri_get_vec(master_table, c("cross_source_robustness"), default = NA_character_)
  coloc_support <- tri_choose_coloc_support(master_table)
  discordance_type <- tri_get_vec(master_table, c("discordance_type"), default = NA_character_)
  technical_caution <- tri_get_vec(master_table, c("technical_caution"), default = NA)

  master_table$triangulation_summary_class <- mapply(
    tri_classify_summary_class,
    eqtl_class = eqtl_class,
    eqtl_primary_support = eqtl_primary_support,
    eqtl_replication_support = eqtl_replication_support,
    pqtl_ukb_support = pqtl_ukb_support,
    pqtl_decode_support = pqtl_decode_support,
    cross_source_robustness = cross_source_robustness,
    coloc_support = coloc_support,
    discordance_type = discordance_type,
    technical_caution = technical_caution,
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  )

  master_table$triangulation_summary_note <- mapply(
    tri_collapse_notes,
    tri_get_vec(master_table, c("cross_source_robustness"), default = NA_character_),
    tri_get_vec(master_table, c("discordance_type"), default = NA_character_),
    tri_get_vec(master_table, c("locus_caution_flag"), default = NA_character_),
    tri_get_vec(master_table, c("coloc_caution_note"), default = NA_character_),
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  )

  master_table
}

tri_build_key_summary <- function(master_table) {
  if (is.null(master_table) || !is.data.frame(master_table) || nrow(master_table) == 0) {
    return(data.frame())
  }

  mt <- master_table
  if (!"triangulation_summary_class" %in% names(mt)) {
    mt <- tri_add_summary_class(mt)
  }

  keep <- intersect(
    c(
      "gene",
      "class.y", "class",
      "triangulation_summary_class",
      "triangulation_summary_note",
      "eqtl_primary_beta",
      "eqtl_primary_p",
      "pqtl_ukbppp_support",
      "pqtl_ukbppp_primary_support",
      "pqtl_ukbppp_primary_beta",
      "pqtl_ukbppp_primary_p",
      "pqtl_decode_support",
      "pqtl_decode_primary_support",
      "pqtl_decode_primary_beta",
      "pqtl_decode_primary_p",
      "ukbppp_decode_concordance",
      "cross_source_robustness",
      "cross_proxy_consistency",
      "discordance_type",
      "final_coloc_support", "coloc_support",
      "locus_caution_flag",
      "manual_review_required",
      "technical_caution",
      "primary_reporting_unit",
      "primary_reporting_label",
      "eqtl_primary_overlap_status",
      "ukbppp_primary_overlap_status",
      "decode_primary_overlap_status"
    ),
    names(mt)
  )

  out <- mt[, keep, drop = FALSE]

  if ("class.y" %in% names(out) && !("eqtl_class" %in% names(out))) {
    names(out)[names(out) == "class.y"] <- "eqtl_class"
  } else if ("class" %in% names(out) && !("eqtl_class" %in% names(out))) {
    names(out)[names(out) == "class"] <- "eqtl_class"
  }

  if ("final_coloc_support" %in% names(out) && !("coloc_support_final" %in% names(out))) {
    names(out)[names(out) == "final_coloc_support"] <- "coloc_support_final"
  } else if ("coloc_support" %in% names(out) && !("coloc_support_final" %in% names(out))) {
    names(out)[names(out) == "coloc_support"] <- "coloc_support_final"
  }

  out
}

tri_build_source_pair_table <- function(pqtl_preferred_df) {
  if (is.null(pqtl_preferred_df) || !is.data.frame(pqtl_preferred_df) || nrow(pqtl_preferred_df) == 0) {
    return(data.frame())
  }

  source_col <- tri_pick_col(pqtl_preferred_df, c("pqtl_source_id", "pqtl_source"))
  if (is.na(source_col)) return(data.frame())

  p1 <- pqtl_preferred_df %>%
    dplyr::filter(.data$outcome_name == "RA_EUR_PRIMARY") %>%
    dplyr::mutate(source_id = .data[[source_col]]) %>%
    dplyr::select(gene, source_id, beta = b, se, pval, method)

  ukb <- p1 %>% dplyr::filter(source_id == "UKB_PPP_EUR") %>%
    dplyr::rename(
      ukb_beta = beta,
      ukb_se = se,
      ukb_p = pval,
      ukb_method = method
    ) %>%
    dplyr::select(-source_id)

  dec <- p1 %>% dplyr::filter(source_id == "deCODE_EUR") %>%
    dplyr::rename(
      decode_beta = beta,
      decode_se = se,
      decode_p = pval,
      decode_method = method
    ) %>%
    dplyr::select(-source_id)

  out <- dplyr::full_join(ukb, dec, by = "gene")
  if (nrow(out) == 0) return(out)

  out$ukb_decode_direction_relation <- dplyr::case_when(
    is.na(out$ukb_beta) | is.na(out$decode_beta) ~ "not_assessable",
    sign(out$ukb_beta) == sign(out$decode_beta) ~ "same_direction",
    sign(out$ukb_beta) != sign(out$decode_beta) ~ "opposite_direction",
    TRUE ~ "not_assessable"
  )

  out
}

tri_build_overlap_snapshot <- function(overlap_dt) {
  if (is.null(overlap_dt) || !is.data.frame(overlap_dt) || nrow(overlap_dt) == 0) {
    return(data.frame())
  }

  keep <- intersect(
    c("exposure_id", "exposure_name", "outcome_name", "outcome_id", "overlap_status", "caveat_level", "recommended_action"),
    names(overlap_dt)
  )

  overlap_dt[, ..keep, drop = FALSE] %>%
    dplyr::distinct() %>%
    dplyr::arrange(.data$exposure_id, .data$outcome_name)
}

tri_write_outputs <- function(dirs,
                              master_table,
                              pqtl_preferred_df = NULL,
                              overlap_dt = NULL,
                              write_master_augmented = FALSE) {
  if (is.null(dirs) || is.null(dirs$results_dir)) {
    stop("dirs$results_dir is required", call. = FALSE)
  }

  mt <- tri_add_summary_class(master_table)

  if (isTRUE(write_master_augmented)) {
    write_csv_safe(mt, fs::path(dirs$results_dir, "master_benchmark_table.csv"))
  }

  write_csv_safe(
    tri_build_key_summary(mt),
    fs::path(dirs$results_dir, "triangulation_key_summary.csv")
  )

  write_csv_safe(
    tri_build_source_pair_table(pqtl_preferred_df),
    fs::path(dirs$results_dir, "triangulation_source_pair_table.csv")
  )

  write_csv_safe(
    tri_build_overlap_snapshot(overlap_dt),
    fs::path(dirs$results_dir, "triangulation_overlap_snapshot.csv")
  )

  invisible(mt)
}

tri_coloc_cross_outcome_concordance <- function(
    coloc_df,
    primary_outcome = "RA_EUR_PRIMARY",
    rep_outcome     = "RA_EUR_REP1_FinnGen") {

  if (is.null(coloc_df) || nrow(coloc_df) == 0) return(data.frame())

  support_col <- tri_pick_col(coloc_df, c("final_coloc_support", "coloc_support"))
  h4_col      <- tri_pick_col(coloc_df, c("final_PP.H4", "PP.H4"))
  gene_col    <- tri_choose_gene_col(coloc_df)

  if (is.na(support_col) || is.na(gene_col)) return(data.frame())

  primary <- coloc_df %>%
    dplyr::filter(outcome_name == primary_outcome) %>%
    dplyr::select(
      gene                  = dplyr::all_of(gene_col),
      primary_coloc_support = dplyr::all_of(support_col),
      primary_pp_h4         = dplyr::all_of(h4_col)
    )

  rep <- coloc_df %>%
    dplyr::filter(outcome_name == rep_outcome) %>%
    dplyr::select(
      gene              = dplyr::all_of(gene_col),
      rep_coloc_support = dplyr::all_of(support_col),
      rep_pp_h4         = dplyr::all_of(h4_col)
    )

  dplyr::full_join(primary, rep, by = "gene") %>%
    dplyr::mutate(
      coloc_cross_outcome_concordance = dplyr::case_when(
        primary_coloc_support == "Strong_coloc" &
          rep_coloc_support == "Strong_coloc"      ~ "Strong_both",
        primary_coloc_support == "Strong_coloc" |
          rep_coloc_support == "Strong_coloc"       ~ "Strong_one_only",
        primary_coloc_support == "Weak_or_uncertain" &
          rep_coloc_support == "Weak_or_uncertain"  ~ "Weak_both",
        is.na(primary_coloc_support) |
          is.na(rep_coloc_support)                  ~ "Not_run_one_or_both",
        TRUE                                         ~ "No_support_both"
      )
    )
}
