# R/11_reporting_tables.R
# Reporting table utilities
# Version: 1.0
# Author: Mihye Kwon
suppressPackageStartupMessages({
  library(dplyr)
})

collapse_notes <- function(...) {
  parts <- unlist(list(...), use.names = FALSE)
  parts <- trimws(as.character(parts))
  parts <- parts[nzchar(parts) & !is.na(parts)]
  if (length(parts) == 0) return(NA_character_)
  paste(unique(parts), collapse = " | ")
}

collect_protocol_note_for_gene <- function(gene, protocol_cfg) {
  notes <- character(0)

  if (!is.null(protocol_cfg$cluster_groups) && length(protocol_cfg$cluster_groups) > 0) {
    for (grp in protocol_cfg$cluster_groups) {
      if (tolower(gene) %in% tolower(grp$genes %||% character(0))) {
        notes <- c(notes, grp$interpretation_rule %||% NA_character_)
      }
    }
  }

  if (!is.null(protocol_cfg$locus_specific_rules) && length(protocol_cfg$locus_specific_rules) > 0) {
    for (rule in protocol_cfg$locus_specific_rules) {
      if (tolower(gene) %in% tolower(rule$genes %||% character(0))) {
        notes <- c(notes, rule$reporting_note %||% rule$interpretation_rule %||% NA_character_)
      }
    }
  }

  collapse_notes(notes)
}



safe_max_numeric <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

gene_locus_caution_flag <- function(gene) {
  g <- tolower(gene)
  if (g %in% c("tnf", "lta")) return("mhc_shared_locus")
  if (g %in% c("fcgr1a", "fcgr2a", "fcgr2b", "fcgr2c", "fcgr3a", "fcgr3b")) return("fcgr_cluster")
  "none"
}

gene_manual_review_required <- function(gene) {
  gene_locus_caution_flag(gene) != "none"
}

gene_technical_caution <- function(gene) {
  gene_manual_review_required(gene)
}

gene_primary_reporting_unit <- function(gene) {
  g <- tolower(gene)
  if (g %in% c("fcgr1a", "fcgr2a", "fcgr2b", "fcgr2c", "fcgr3a", "fcgr3b")) {
    return("cluster_level_primary")
  }
  if (g %in% c("tnf", "lta")) {
    return("gene_level_with_shared_locus_caution")
  }
  "gene_level"
}

gene_primary_reporting_label <- function(gene) {
  g <- tolower(gene)
  if (g %in% c("fcgr1a", "fcgr2a", "fcgr2b", "fcgr2c", "fcgr3a", "fcgr3b")) {
    return("FCGR_cluster_1q23")
  }
  if (g %in% c("tnf", "lta")) {
    return("TNF_LTA_shared_locus")
  }
  toupper(g)
}

build_coloc_caution_note <- function(gene,
                                     coloc_support = NA_character_,
                                     final_coloc_support = NA_character_,
                                     susie_status = NA_character_,
                                     ld_matrix_note = NA_character_,
                                     note = NA_character_) {
  g <- tolower(gene)
  out <- character(0)

  final_support <- collapse_notes(final_coloc_support, coloc_support)
  if (!is.na(final_support) && final_support %in% c("Distinct_signals_likely", "Weak_or_uncertain")) {
    out <- c(out, final_support)
  }
  if (!is.na(susie_status) && susie_status == "failed") {
    out <- c(out, "coloc.susie_failed")
  }
  if (g %in% c("tnf", "lta")) {
    out <- c(out, "Complex shared locus / MHC-proximal caution")
  }
  collapse_notes(out, ld_matrix_note, note)
}

classify_pqtl_estimability_state <- function(df) {
  if (is.null(df) || nrow(df) == 0) return("not_applicable")
  if (any(isTRUE(df$estimable), na.rm = TRUE)) return("estimable")

  reasons <- unique(stats::na.omit(as.character(df$non_estimable_reason)))
  raw_counts <- suppressWarnings(as.numeric(df$raw_variant_count))

  if (length(raw_counts) > 0 && all(is.na(raw_counts) | raw_counts == 0) &&
      all(reasons %in% c("no_variant_after_target_filter"))) {
    return("source_unavailable")
  }

  "non_estimable"
}

summarise_single_pqtl_source_attrition <- function(pqtl_attrition_df,
                                                   target_genes,
                                                   source_id,
                                                   prefix) {
  base_tbl <- data.frame(gene = tolower(target_genes), stringsAsFactors = FALSE)

  if (is.null(pqtl_attrition_df) || nrow(pqtl_attrition_df) == 0 ||
      !source_id %in% pqtl_attrition_df$pqtl_source_id) {
    base_tbl[[paste0(prefix, "_estimability")]] <- "not_applicable"
    base_tbl[[paste0(prefix, "_non_estimable_reason")]] <- NA_character_
    base_tbl[[paste0(prefix, "_raw_variant_count")]] <- NA_real_
    base_tbl[[paste0(prefix, "_post_p_threshold_count")]] <- NA_real_
    base_tbl[[paste0(prefix, "_post_clump_count")]] <- NA_real_
    base_tbl[[paste0(prefix, "_post_harmonisation_usable_count")]] <- NA_real_
    return(base_tbl)
  }

  dat <- pqtl_attrition_df %>% filter(pqtl_source_id == source_id)
  if ("outcome_name" %in% names(dat) && any(dat$outcome_name == "RA_EUR_PRIMARY")) {
    dat <- dat %>% filter(outcome_name == "RA_EUR_PRIMARY")
  }

  sumdat <- dat %>%
    group_by(gene) %>%
    summarise(
      estimability = classify_pqtl_estimability_state(cur_data_all()),
      non_estimable_reason = collapse_notes(non_estimable_reason),
      raw_variant_count = safe_max_numeric(raw_variant_count),
      post_p_threshold_count = safe_max_numeric(post_p_threshold_count),
      post_clump_count = safe_max_numeric(post_clump_count),
      post_harmonisation_usable_count = safe_max_numeric(post_harmonisation_usable_count),
      .groups = "drop"
    )

  names(sumdat)[names(sumdat) == "estimability"] <- paste0(prefix, "_estimability")
  names(sumdat)[names(sumdat) == "non_estimable_reason"] <- paste0(prefix, "_non_estimable_reason")
  names(sumdat)[names(sumdat) == "raw_variant_count"] <- paste0(prefix, "_raw_variant_count")
  names(sumdat)[names(sumdat) == "post_p_threshold_count"] <- paste0(prefix, "_post_p_threshold_count")
  names(sumdat)[names(sumdat) == "post_clump_count"] <- paste0(prefix, "_post_clump_count")
  names(sumdat)[names(sumdat) == "post_harmonisation_usable_count"] <- paste0(prefix, "_post_harmonisation_usable_count")

  base_tbl %>% left_join(sumdat, by = "gene") %>%
    mutate(
      !!paste0(prefix, "_estimability") := .data[[paste0(prefix, "_estimability")]] %||% "not_applicable"
    )
}


summarise_single_pqtl_source_support <- function(preferred_pqtl_df,
                                                 target_meta,
                                                 target_genes,
                                                 source_id,
                                                 prefix,
                                                 nominal_p_cut = 0.05) {
  base_tbl <- data.frame(gene = tolower(target_genes), stringsAsFactors = FALSE)

  base_tbl[[paste0(prefix, "_support")]] <- "Not_estimable"
  base_tbl[[paste0(prefix, "_primary_support")]] <- FALSE
  base_tbl[[paste0(prefix, "_replication_support")]] <- FALSE
  base_tbl[[paste0(prefix, "_primary_beta")]] <- NA_real_
  base_tbl[[paste0(prefix, "_primary_p")]] <- NA_real_
  base_tbl[[paste0(prefix, "_n_rows")]] <- 0L

  if (is.null(preferred_pqtl_df) || nrow(preferred_pqtl_df) == 0) {
    return(base_tbl)
  }

  source_col <- if ("pqtl_source_id" %in% names(preferred_pqtl_df)) {
    "pqtl_source_id"
  } else if ("pqtl_source" %in% names(preferred_pqtl_df)) {
    "pqtl_source"
  } else {
    NULL
  }

  if (is.null(source_col)) {
    return(base_tbl)
  }

  src_dat <- preferred_pqtl_df %>% dplyr::filter(.data[[source_col]] == source_id)

  if (nrow(src_dat) == 0) {
    return(base_tbl)
  }

  out <- lapply(tolower(target_genes), function(g) {
    gdat <- src_dat %>% dplyr::filter(.data$gene == g)
    tm <- if (!is.null(target_meta) && nrow(target_meta) > 0) {
      target_meta %>% dplyr::filter(.data$gene == g)
    } else {
      NULL
    }

    p1 <- gdat %>% dplyr::filter(.data$outcome_name == "RA_EUR_PRIMARY")
    rp <- if ("role" %in% names(gdat)) {
      gdat %>% dplyr::filter(.data$role == "replication")
    } else {
      gdat[0, , drop = FALSE]
    }

    primary_support_vec <- logical(0)
    primary_sign_vec <- numeric(0)

    if (nrow(p1) > 0) {
      primary_support_vec <- vapply(seq_len(nrow(p1)), function(i) {
        dir_mod <- expected_direction_modifier(tm, p1$b[i])
        isTRUE(safe_sig(p1[i, , drop = FALSE], nominal_p_cut)) &&
          !identical(dir_mod$expected_direction_concordant, FALSE)
      }, logical(1))
      primary_sign_vec <- p1$b
    }

    replication_support <- FALSE
    if (length(primary_sign_vec) > 0 && nrow(rp) > 0) {
      replication_support <- any(vapply(seq_len(nrow(rp)), function(i) {
        isTRUE(safe_sig(rp[i, , drop = FALSE], nominal_p_cut)) &&
          any(vapply(primary_sign_vec, function(pb) same_dir(pb, rp$b[i]), logical(1)), na.rm = TRUE)
      }, logical(1)), na.rm = TRUE)
    }

    primary_support_any <- if (length(primary_support_vec) > 0) any(primary_support_vec, na.rm = TRUE) else FALSE

    support_label <- if (nrow(gdat) == 0) {
      "Not_estimable"
    } else if (isTRUE(primary_support_any) && isTRUE(replication_support)) {
      "Supportive_primary_and_replication"
    } else if (isTRUE(primary_support_any)) {
      "Supportive_primary_only"
    } else {
      "No_support_detected"
    }

    data.frame(
      gene = g,
      support = support_label,
      primary_support = primary_support_any,
      replication_support = replication_support,
      primary_beta = if (nrow(p1) > 0) p1$b[1] else NA_real_,
      primary_p = if (nrow(p1) > 0) p1$pval[1] else NA_real_,
      n_rows = nrow(gdat),
      stringsAsFactors = FALSE
    )
  })

  out <- bind_rows(out)
  names(out)[names(out) == "support"] <- paste0(prefix, "_support")
  names(out)[names(out) == "primary_support"] <- paste0(prefix, "_primary_support")
  names(out)[names(out) == "replication_support"] <- paste0(prefix, "_replication_support")
  names(out)[names(out) == "primary_beta"] <- paste0(prefix, "_primary_beta")
  names(out)[names(out) == "primary_p"] <- paste0(prefix, "_primary_p")
  names(out)[names(out) == "n_rows"] <- paste0(prefix, "_n_rows")

  base_tbl %>%
    dplyr::select(gene) %>%
    dplyr::left_join(out, by = "gene")
}

build_overlap_wide <- function(overlap_dt) {
  if (is.null(overlap_dt) || !is.data.frame(overlap_dt) || nrow(overlap_dt) == 0) {
    return(data.frame())
  }

  required_cols <- c("exposure_id", "outcome_name", "overlap_status")
  if (!all(required_cols %in% names(overlap_dt))) {
    return(data.frame())
  }

  overlap_map <- data.frame(
    exposure_id = c(
      "EQTLGEN_BLOOD", "EQTLGEN_BLOOD", "EQTLGEN_BLOOD", "EQTLGEN_BLOOD", "EQTLGEN_BLOOD",
      "UKB_PPP_EUR", "UKB_PPP_EUR", "UKB_PPP_EUR", "UKB_PPP_EUR", "UKB_PPP_EUR",
      "deCODE_EUR", "deCODE_EUR", "deCODE_EUR", "deCODE_EUR", "deCODE_EUR"
    ),
    outcome_name = c(
      "RA_EUR_PRIMARY", "RA_EUR_REP1_FinnGen", "RA_EUR_SEROPOS", "RA_EUR_SERONEG", "RA_EAS_BBJ",
      "RA_EUR_PRIMARY", "RA_EUR_REP1_FinnGen", "RA_EUR_SEROPOS", "RA_EUR_SERONEG", "RA_EAS_BBJ",
      "RA_EUR_PRIMARY", "RA_EUR_REP1_FinnGen", "RA_EUR_SEROPOS", "RA_EUR_SERONEG", "RA_EAS_BBJ"
    ),
    out_col = c(
      "eqtl_primary_overlap_status", "eqtl_replication_overlap_status", "eqtl_seropos_overlap_status", "eqtl_seroneg_overlap_status", "eqtl_eas_overlap_status",
      "ukbppp_primary_overlap_status", "ukbppp_replication_overlap_status", "ukbppp_seropos_overlap_status", "ukbppp_seroneg_overlap_status", "ukbppp_eas_overlap_status",
      "decode_primary_overlap_status", "decode_replication_overlap_status", "decode_seropos_overlap_status", "decode_seroneg_overlap_status", "decode_eas_overlap_status"
    ),
    stringsAsFactors = FALSE
  )

  ov <- overlap_dt %>%
    dplyr::select(exposure_id, outcome_name, overlap_status) %>%
    dplyr::distinct() %>%
    dplyr::left_join(overlap_map, by = c("exposure_id", "outcome_name")) %>%
    dplyr::filter(!is.na(.data$out_col))

  out <- as.list(rep(NA_character_, length(unique(overlap_map$out_col))))
  names(out) <- unique(overlap_map$out_col)

  if (nrow(ov) > 0) {
    for (i in seq_len(nrow(ov))) {
      out[[ov$out_col[i]]] <- ov$overlap_status[i]
    }
  }

  as.data.frame(out, stringsAsFactors = FALSE)
}

classify_cross_source_robustness <- function(ukb_state,
                                             decode_state,
                                             ukb_primary_support,
                                             decode_primary_support,
                                             ukb_beta,
                                             decode_beta) {
  if (isTRUE(ukb_primary_support) && isTRUE(decode_primary_support)) {
    if (!is.na(ukb_beta) && !is.na(decode_beta)) {
      if (sign(ukb_beta) == sign(decode_beta)) {
        return("concordant_support_two_source")
      } else {
        return("discordant_support_two_source")
      }
    }
    return("two_source_support_direction_unknown")
  }

  if (xor(isTRUE(ukb_primary_support), isTRUE(decode_primary_support))) {
    return("single_source_support")
  }

  if (identical(ukb_state, "estimable") && identical(decode_state, "estimable")) {
    return("two_source_estimable_no_dual_support")
  }

  if (xor(identical(ukb_state, "estimable"), identical(decode_state, "estimable"))) {
    return("single_source_only")
  }

  if (identical(ukb_state, "source_unavailable") && identical(decode_state, "source_unavailable")) {
    return("source_unavailable_both")
  }

  "not_assessable"
}

summarise_pqtl_support <- function(preferred_pqtl_df, target_meta = NULL, nominal_p_cut = 0.05) {
  if (is.null(preferred_pqtl_df) || nrow(preferred_pqtl_df) == 0) return(data.frame())

  genes <- sort(unique(preferred_pqtl_df$gene))

  out <- lapply(genes, function(g) {
    gdat <- preferred_pqtl_df %>% filter(gene == g)
    p1 <- gdat %>% filter(outcome_name == "RA_EUR_PRIMARY")
    rp <- gdat %>% filter(role == "replication")
    tm <- if (!is.null(target_meta)) target_meta %>% filter(gene == g) else NULL

    if (nrow(gdat) == 0) {
      return(data.frame(
        gene = g,
        pqtl_support = "Not_estimable",
        pqtl_estimable_any = FALSE,
        pqtl_primary_support_any = FALSE,
        pqtl_replication_support_any = FALSE,
        pqtl_n_sources = 0L,
        pqtl_n_sources_primary_support = 0L,
        stringsAsFactors = FALSE
      ))
    }

    primary_support_vec <- logical(0)
    primary_sign_vec <- numeric(0)

    if (nrow(p1) > 0) {
      primary_support_vec <- vapply(seq_len(nrow(p1)), function(i) {
        dir_mod <- expected_direction_modifier(tm, p1$b[i])
        isTRUE(safe_sig(p1[i, , drop = FALSE], nominal_p_cut)) &&
          !identical(dir_mod$expected_direction_concordant, FALSE)
      }, logical(1))
      primary_sign_vec <- p1$b
    }

    replication_support <- FALSE
    if (length(primary_sign_vec) > 0 && nrow(rp) > 0) {
      replication_support <- any(vapply(seq_len(nrow(rp)), function(i) {
        isTRUE(safe_sig(rp[i, , drop = FALSE], nominal_p_cut)) &&
          any(vapply(primary_sign_vec, function(pb) same_dir(pb, rp$b[i]), logical(1)), na.rm = TRUE)
      }, logical(1)), na.rm = TRUE)
    }

    primary_support_any <- if (length(primary_support_vec) > 0) any(primary_support_vec, na.rm = TRUE) else FALSE
    estimable_any <- nrow(gdat) > 0
    pqtl_support <- if (!estimable_any) {
      "Not_estimable"
    } else if (isTRUE(primary_support_any) && isTRUE(replication_support)) {
      "Supportive_primary_and_replication"
    } else if (isTRUE(primary_support_any)) {
      "Supportive_primary_only"
    } else {
      "No_support_detected"
    }

    data.frame(
      gene = g,
      pqtl_support = pqtl_support,
      pqtl_estimable_any = estimable_any,
      pqtl_primary_support_any = primary_support_any,
      pqtl_replication_support_any = replication_support,
      pqtl_n_sources = dplyr::n_distinct(gdat$pqtl_source),
      pqtl_n_sources_primary_support = if (nrow(p1) > 0) sum(primary_support_vec, na.rm = TRUE) else 0L,
      stringsAsFactors = FALSE
    )
  })

  bind_rows(out)
}

build_cross_proxy_consistency_table <- function(eqtl_preferred_df, pqtl_preferred_df) {
  if (is.null(eqtl_preferred_df) || nrow(eqtl_preferred_df) == 0) return(data.frame())

  eqtl_primary <- eqtl_preferred_df %>%
    filter(outcome_name == "RA_EUR_PRIMARY") %>%
    select(gene, eqtl_beta = b, eqtl_p = pval) %>%
    mutate(
      eqtl_primary_available = !is.na(eqtl_beta),
      eqtl_primary_significant = !is.na(eqtl_p) & eqtl_p < 0.05
    )

  if (is.null(pqtl_preferred_df) || nrow(pqtl_preferred_df) == 0) {
    return(eqtl_primary %>%
      mutate(
        pqtl_primary_available = FALSE,
        pqtl_sources_available = 0L,
        pqtl_sources_same_direction = 0L,
        pqtl_sources_opposite_direction = 0L,
        cross_proxy_assessable = FALSE,
        cross_proxy_consistency = "No_pQTL_data",
        discordance_type = "Type1_no_pQTL_data"
      ))
  }

  pqtl_primary <- pqtl_preferred_df %>%
    filter(outcome_name == "RA_EUR_PRIMARY") %>%
    select(gene, pqtl_source, pqtl_beta = b, pqtl_p = pval)

  full_join(eqtl_primary, pqtl_primary, by = "gene") %>%
    group_by(gene, eqtl_beta, eqtl_p, eqtl_primary_available, eqtl_primary_significant) %>%
    summarise(
      pqtl_primary_available = any(!is.na(pqtl_beta)),
      pqtl_sources_available = sum(!is.na(pqtl_beta)),
      pqtl_sources_same_direction = sum(sign(eqtl_beta) == sign(pqtl_beta), na.rm = TRUE),
      pqtl_sources_opposite_direction = sum(sign(eqtl_beta) != sign(pqtl_beta), na.rm = TRUE),
      cross_proxy_assessable = isTRUE(first(eqtl_primary_available)) && any(!is.na(pqtl_beta)),
      cross_proxy_consistency = dplyr::case_when(
        !any(!is.na(pqtl_beta)) ~ "No_pQTL_data",
        sum(sign(eqtl_beta) == sign(pqtl_beta), na.rm = TRUE) > 0 &&
          sum(sign(eqtl_beta) != sign(pqtl_beta), na.rm = TRUE) == 0 ~ "Concordant",
        sum(sign(eqtl_beta) == sign(pqtl_beta), na.rm = TRUE) == 0 &&
          sum(sign(eqtl_beta) != sign(pqtl_beta), na.rm = TRUE) > 0 ~ "Discordant",
        TRUE ~ "Mixed"
      ),
      discordance_type = dplyr::case_when(
        !any(!is.na(pqtl_beta)) ~ "Type1_no_pQTL_data",
        !isTRUE(first(eqtl_primary_significant)) ~ "Type4_eqtl_null_pqtl_present",
        sum(sign(eqtl_beta) == sign(pqtl_beta), na.rm = TRUE) == 0 &&
          sum(sign(eqtl_beta) != sign(pqtl_beta), na.rm = TRUE) > 0 ~ "Type2_eqtl_pqtl_opposite",
        TRUE ~ NA_character_
      ),
      .groups = "drop"
    )
}

augment_classification_support_layers <- function(class_df,
                                                  target_meta,
                                                  protocol_cfg,
                                                  coloc_df = NULL,
                                                  pqtl_support_df = NULL,
                                                  cross_proxy_df = NULL) {
  out <- class_df

  if (!is.null(coloc_df) && nrow(coloc_df) > 0) {
    coloc_join <- coloc_df
    if ("outcome_name" %in% names(coloc_join) &&
        any(coloc_join$outcome_name == "RA_EUR_PRIMARY")) {
      coloc_join <- coloc_join %>%
        dplyr::filter(outcome_name == "RA_EUR_PRIMARY")
    }
    coloc_join <- coloc_join %>%
      mutate(
        coloc_support_join = if ("final_coloc_support" %in% names(.)) .data$final_coloc_support else .data$coloc_support,
        coloc_pp_h4_join = if ("final_PP.H4" %in% names(.)) .data$final_PP.H4 else .data$PP.H4
      ) %>%
      select(gene, coloc_support = coloc_support_join, coloc_pp_h4 = coloc_pp_h4_join)
    out <- out %>% left_join(coloc_join, by = "gene")
  } else {
    out$coloc_support <- NA_character_
    out$coloc_pp_h4 <- NA_real_
  }

  if (!is.null(pqtl_support_df) && nrow(pqtl_support_df) > 0) {
    out <- out %>% left_join(pqtl_support_df, by = "gene")
  } else {
    out$pqtl_support <- NA_character_
    out$pqtl_estimable_any <- NA
    out$pqtl_primary_support_any <- NA
    out$pqtl_replication_support_any <- NA
    out$pqtl_n_sources <- NA_integer_
    out$pqtl_n_sources_primary_support <- NA_integer_
  }

  if (!is.null(cross_proxy_df) && nrow(cross_proxy_df) > 0) {
    out <- out %>% left_join(
      cross_proxy_df %>% select(
        gene, cross_proxy_consistency, pqtl_primary_available,
        pqtl_sources_available, pqtl_sources_same_direction, pqtl_sources_opposite_direction,
        discordance_type
      ),
      by = "gene"
    )
  } else {
    out$cross_proxy_consistency <- NA_character_
    out$pqtl_primary_available <- NA
    out$pqtl_sources_available <- NA_integer_
    out$pqtl_sources_same_direction <- NA_integer_
    out$pqtl_sources_opposite_direction <- NA_integer_
    out$discordance_type <- NA_character_
  }

  if (!is.null(target_meta) && nrow(target_meta) > 0) {
    out <- out %>%
      left_join(
        target_meta %>%
          select(
            gene,
            tier,
            target_role,
            target_class = class,
            drug_examples,
            interpretability_note,
            alignment_certainty
          ),
        by = "gene"
      )
  } else {
    out$tier <- NA_character_
    out$target_role <- NA_character_
    out$target_class <- NA_character_
    out$drug_examples <- NA_character_
    out$interpretability_note <- NA_character_
    out$alignment_certainty <- NA_character_
  }

  out$protocol_note <- vapply(
    out$gene,
    collect_protocol_note_for_gene,
    character(1),
    protocol_cfg = protocol_cfg
  )

  out$interpretation_note <- mapply(
    collapse_notes,
    out$interpretability_note %||% NA_character_,
    out$protocol_note %||% NA_character_,
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  )

  out
}

build_master_benchmark_table <- function(target_meta,
                                         attrition_df,
                                         eqtl_preferred_df,
                                         sero_df,
                                         transport_df,
                                         class_support_df,
                                         coloc_df = NULL,
                                         pqtl_support_df = NULL,
                                         cross_proxy_df = NULL,
                                         pqtl_attrition_df = NULL,
                                         pqtl_preferred_df = NULL,
                                         overlap_dt = NULL) {
  target_tbl <- target_meta %>%
    dplyr::select(
      gene, tier, target_role, class, direct_ra_target,
      drug_examples, expected_protective_direction_if_aligned,
      alignment_certainty, interpretability_note
    )

  if (!"nsnp" %in% names(eqtl_preferred_df) && "nsnp_usable" %in% names(eqtl_preferred_df)) {
    eqtl_preferred_df <- eqtl_preferred_df %>% dplyr::rename(nsnp = nsnp_usable)
  }

  primary_eqtl <- eqtl_preferred_df %>%
    dplyr::filter(outcome_name == "RA_EUR_PRIMARY") %>%
    dplyr::select(
      gene,
      eqtl_primary_method = method,
      eqtl_primary_beta = b,
      eqtl_primary_se = se,
      eqtl_primary_p = pval,
      eqtl_primary_nsnp = nsnp
    )

  rep_eqtl <- eqtl_preferred_df %>%
    dplyr::filter(role == "replication") %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      eqtl_replication_any = n() > 0,
      eqtl_replication_beta_first = dplyr::first(b),
      eqtl_replication_p_first = dplyr::first(pval),
      .groups = "drop"
    )

  primary_attr <- attrition_df %>%
    dplyr::filter(outcome_name == "RA_EUR_PRIMARY") %>%
    dplyr::select(
      gene, raw_variant_count, post_p_threshold_count, post_clump_count,
      post_harmonisation_usable_count, estimable, non_estimable_reason
    )

  out <- target_tbl %>%
    dplyr::left_join(primary_attr, by = "gene") %>%
    dplyr::left_join(primary_eqtl, by = "gene") %>%
    dplyr::left_join(rep_eqtl, by = "gene") %>%
    dplyr::left_join(sero_df, by = "gene") %>%
    dplyr::left_join(transport_df, by = "gene") %>%
    dplyr::left_join(class_support_df, by = "gene")

  if (!is.null(coloc_df) && nrow(coloc_df) > 0) {
    coloc_join <- coloc_df
    if ("outcome_name" %in% names(coloc_join) && any(coloc_join$outcome_name == "RA_EUR_PRIMARY")) {
      coloc_join <- coloc_join %>% dplyr::filter(outcome_name == "RA_EUR_PRIMARY")
    }

    keep_cols <- intersect(
      c("gene", "PP.H4", "PP.H3", "coloc_support", "final_coloc_support",
        "susie_status", "ld_matrix_note", "note"),
      names(coloc_join)
    )
    coloc_join <- coloc_join[, keep_cols, drop = FALSE]

    if (!"final_coloc_support" %in% names(coloc_join)) coloc_join$final_coloc_support <- NA_character_
    if (!"susie_status" %in% names(coloc_join)) coloc_join$susie_status <- NA_character_
    if (!"ld_matrix_note" %in% names(coloc_join)) coloc_join$ld_matrix_note <- NA_character_
    if (!"note" %in% names(coloc_join)) coloc_join$note <- NA_character_

    coloc_join$coloc_caution_note <- mapply(
      build_coloc_caution_note,
      gene = coloc_join$gene,
      coloc_support = coloc_join$coloc_support %||% NA_character_,
      final_coloc_support = coloc_join$final_coloc_support %||% NA_character_,
      susie_status = coloc_join$susie_status %||% NA_character_,
      ld_matrix_note = coloc_join$ld_matrix_note %||% NA_character_,
      note = coloc_join$note %||% NA_character_,
      SIMPLIFY = TRUE,
      USE.NAMES = FALSE
    )

    out <- out %>% dplyr::left_join(coloc_join, by = "gene")
  }

  if (!is.null(pqtl_support_df) && nrow(pqtl_support_df) > 0) {
    out <- out %>% dplyr::left_join(pqtl_support_df, by = "gene", suffix = c("", ".pqtl"))
  }

  if (!is.null(cross_proxy_df) && nrow(cross_proxy_df) > 0) {
    out <- out %>% dplyr::left_join(cross_proxy_df, by = "gene", suffix = c("", ".cross"))
  }

  ukbppp_state <- summarise_single_pqtl_source_attrition(
    pqtl_attrition_df = pqtl_attrition_df,
    target_genes = target_tbl$gene,
    source_id = "UKB_PPP_EUR",
    prefix = "pqtl_ukbppp"
  )

  decode_state <- summarise_single_pqtl_source_attrition(
    pqtl_attrition_df = pqtl_attrition_df,
    target_genes = target_tbl$gene,
    source_id = "deCODE_EUR",
    prefix = "pqtl_decode"
  )

  ukbppp_support <- summarise_single_pqtl_source_support(
    preferred_pqtl_df = pqtl_preferred_df,
    target_meta = target_meta,
    target_genes = target_tbl$gene,
    source_id = "UKB_PPP_EUR",
    prefix = "pqtl_ukbppp"
  )

  decode_support <- summarise_single_pqtl_source_support(
    preferred_pqtl_df = pqtl_preferred_df,
    target_meta = target_meta,
    target_genes = target_tbl$gene,
    source_id = "deCODE_EUR",
    prefix = "pqtl_decode"
  )

  out <- out %>%
    dplyr::left_join(ukbppp_state, by = "gene") %>%
    dplyr::left_join(decode_state, by = "gene") %>%
    dplyr::left_join(ukbppp_support, by = "gene") %>%
    dplyr::left_join(decode_support, by = "gene")

  out$locus_caution_flag <- vapply(out$gene, gene_locus_caution_flag, character(1))
  out$manual_review_required <- vapply(out$gene, gene_manual_review_required, logical(1))
  out$technical_caution <- vapply(out$gene, gene_technical_caution, logical(1))
  out$primary_reporting_unit <- vapply(out$gene, gene_primary_reporting_unit, character(1))
  out$primary_reporting_label <- vapply(out$gene, gene_primary_reporting_label, character(1))
  out$cluster_primary_reporting_rule <- ifelse(
    out$primary_reporting_unit == "cluster_level_primary",
    "main_text_cluster_level__supplement_gene_level_sensitivity",
    NA_character_
  )

  out$ukbppp_decode_concordance <- dplyr::case_when(
    is.na(out$pqtl_ukbppp_primary_beta) | is.na(out$pqtl_decode_primary_beta) ~ "not_assessable",
    sign(out$pqtl_ukbppp_primary_beta) == sign(out$pqtl_decode_primary_beta) ~ "concordant",
    sign(out$pqtl_ukbppp_primary_beta) != sign(out$pqtl_decode_primary_beta) ~ "discordant",
    TRUE ~ "not_assessable"
  )

  source_type3_discordance <- dplyr::case_when(
    isTRUE(out$pqtl_ukbppp_primary_support) &
      isTRUE(out$pqtl_decode_primary_support) &
      out$ukbppp_decode_concordance == "discordant" ~ "Type3_pqtl_ukbppp_decode_opposite",
    TRUE ~ NA_character_
  )

  out$discordance_type <- mapply(
    collapse_notes,
    out$discordance_type %||% NA_character_,
    source_type3_discordance %||% NA_character_,
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  )

  out$cross_source_robustness <- mapply(
    classify_cross_source_robustness,
    ukb_state = out$pqtl_ukbppp_estimability,
    decode_state = out$pqtl_decode_estimability,
    ukb_primary_support = out$pqtl_ukbppp_primary_support,
    decode_primary_support = out$pqtl_decode_primary_support,
    ukb_beta = out$pqtl_ukbppp_primary_beta,
    decode_beta = out$pqtl_decode_primary_beta,
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  )

  if (!is.null(overlap_dt) && nrow(overlap_dt) > 0) {
    ov_wide <- build_overlap_wide(overlap_dt)
    if (nrow(ov_wide) > 0) {
      out <- dplyr::bind_cols(
        out,
        ov_wide[rep(1, nrow(out)), , drop = FALSE]
      )
    }
  }

  out
}

build_replication_coloc_summary <- function(
    coloc_df,
    replication_outcome_name = "RA_EUR_REP1_FinnGen") {

  if (is.null(coloc_df) || nrow(coloc_df) == 0) return(data.frame())

  keep_cols <- intersect(
    c("gene", "outcome_name", "outcome_id",
      "nsnps_eqtl", "nsnps_overlap",
      "PP.H3", "PP.H4",
      "final_PP.H4", "final_coloc_support",
      "susie_triggered", "susie_status",
      "ld_matrix_note", "note"),
    names(coloc_df)
  )

  coloc_df %>%
    dplyr::filter(outcome_name == replication_outcome_name) %>%
    dplyr::select(dplyr::all_of(keep_cols))
}
