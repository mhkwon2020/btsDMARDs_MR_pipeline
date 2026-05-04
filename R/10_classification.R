# R/10_classification.R
# Benchmark classification utilities
# Version: 1.0
# Author: Mihye Kwon
suppressPackageStartupMessages({
  library(dplyr)
})

safe_sig <- function(df, p_cut = 0.05) {
  nrow(df) > 0 &&
    !is.na(df$pval[1]) &&
    is.finite(df$pval[1]) &&
    !is.na(df$b[1]) &&
    is.finite(df$b[1]) &&
    df$pval[1] < p_cut
}

same_dir <- function(a, b) {
  !is.na(a) && is.finite(a) &&
    !is.na(b) && is.finite(b) &&
    sign(a) == sign(b)
}

build_preferred_results <- function(all_mr_results, run_cfg) {
  if (is.null(all_mr_results) || !is.data.frame(all_mr_results) || nrow(all_mr_results) == 0) {
    return(data.frame(
      gene = character(0),
      outcome_name = character(0),
      outcome_id = character(0),
      role = character(0),
      ancestry = character(0),
      phenotype = character(0),
      method = character(0),
      nsnp_usable = integer(0),
      b = numeric(0),
      se = numeric(0),
      pval = numeric(0),
      OR = numeric(0),
      OR_lci95 = numeric(0),
      OR_uci95 = numeric(0),
      f_median = numeric(0),
      f_min = numeric(0),
      f_max = numeric(0),
      n_f_gt_10 = integer(0),
      stringsAsFactors = FALSE
    ))
  }
  method_order <- preferred_method_order(run_cfg)
  split_df <- split(all_mr_results, list(all_mr_results$gene, all_mr_results$outcome_name), drop = TRUE)
  preferred <- lapply(split_df, function(df) pick_preferred_result(df, method_order))
  dplyr::bind_rows(preferred)
}

make_serostatus_table <- function(preferred_df) {
  empty <- data.frame(
    gene = character(0),
    b_pos = numeric(0),
    p_pos = numeric(0),
    b_neg = numeric(0),
    p_neg = numeric(0),
    dir_concord = logical(0),
    stringsAsFactors = FALSE
  )
  if (is.null(preferred_df) || nrow(preferred_df) == 0) return(empty)
  if (!"outcome_name" %in% names(preferred_df)) return(empty)

  pos <- preferred_df %>% filter(outcome_name == "RA_EUR_SEROPOS") %>%
    select(gene, b_pos = b, p_pos = pval)
  neg <- preferred_df %>% filter(outcome_name == "RA_EUR_SERONEG") %>%
    select(gene, b_neg = b, p_neg = pval)

  if (nrow(pos) == 0 && nrow(neg) == 0) return(empty)

  full_join(pos, neg, by = "gene") %>%
    mutate(
      dir_concord = ifelse(!is.na(b_pos) & !is.na(b_neg), sign(b_pos) == sign(b_neg), NA)
    )
}

make_transportability_table <- function(preferred_df) {
  empty <- data.frame(
    gene = character(0),
    b_eur = numeric(0),
    p_eur = numeric(0),
    b_eas = numeric(0),
    p_eas = numeric(0),
    dir_concord = logical(0),
    stringsAsFactors = FALSE
  )
  if (is.null(preferred_df) || nrow(preferred_df) == 0) return(empty)
  if (!"outcome_name" %in% names(preferred_df)) return(empty)

  eur <- preferred_df %>% filter(outcome_name == "RA_EUR_PRIMARY") %>%
    select(gene, b_eur = b, p_eur = pval)
  eas <- preferred_df %>% filter(outcome_name == "RA_EAS_BBJ") %>%
    select(gene, b_eas = b, p_eas = pval)

  if (nrow(eur) == 0 && nrow(eas) == 0) return(empty)

  full_join(eur, eas, by = "gene") %>%
    mutate(
      dir_concord = ifelse(!is.na(b_eur) & !is.na(b_eas), sign(b_eur) == sign(b_eas), NA)
    )
}

expected_direction_modifier <- function(target_row, beta_value) {
  if (is.null(target_row) || nrow(target_row) == 0 || is.na(beta_value) || !is.finite(beta_value)) {
    return(list(
      expected_direction_concordant = NA,
      alignment_caution_flag = NA_character_
    ))
  }

  exp_dir <- target_row$expected_protective_direction_if_aligned[1]
  align_cert <- target_row$alignment_certainty[1]

  ambiguous_dirs <- c("complex", "context_dependent")
  negative_beta_supportive <- c(
    "lower_target_activity_or_expression",
    "lower_target_availability_or_costimulation",
    "lower_target_cell_abundance_or_target_presence",
    "lower_signalling",
    "higher_inhibitory_signalling"
  )
  positive_beta_supportive <- c(
    "higher_target_activity_or_expression",
    "higher_target_availability_or_costimulation",
    "higher_target_cell_abundance_or_target_presence",
    "higher_signalling"
  )

  if (is.na(exp_dir) || exp_dir %in% ambiguous_dirs) {
    return(list(
      expected_direction_concordant = NA,
      alignment_caution_flag = ifelse(is.na(align_cert), NA_character_, paste0("alignment_", align_cert))
    ))
  }

  concord <- if (exp_dir %in% negative_beta_supportive) {
    beta_value < 0
  } else if (exp_dir %in% positive_beta_supportive) {
    beta_value > 0
  } else {
    NA
  }

  list(
    expected_direction_concordant = concord,
    alignment_caution_flag = ifelse(is.na(align_cert), NA_character_, paste0("alignment_", align_cert))
  )
}

first_or_na <- function(df, col) {
  if (is.null(df) || nrow(df) == 0 || !(col %in% names(df))) return(NA)
  df[[col]][1]
}

make_benchmark_classification <- function(preferred_df, attrition_df, run_cfg, target_meta = NULL) {
  if (is.null(preferred_df) || !is.data.frame(preferred_df)) {
    preferred_df <- data.frame()
  }
  if (!"gene" %in% names(preferred_df))         preferred_df$gene         <- character(0)
  if (!"outcome_name" %in% names(preferred_df)) preferred_df$outcome_name <- character(0)
  if (!"role" %in% names(preferred_df))         preferred_df$role         <- character(0)
  if (!"b" %in% names(preferred_df))            preferred_df$b            <- numeric(0)
  if (!"pval" %in% names(preferred_df))         preferred_df$pval         <- numeric(0)

  genes <- sort(unique(c(attrition_df$gene, preferred_df$gene)))
  nominal_p_cut <- 0.05
  if (!is.null(run_cfg$classification$primary_p_threshold) &&
      is.numeric(run_cfg$classification$primary_p_threshold) &&
      length(run_cfg$classification$primary_p_threshold) == 1 &&
      !is.na(run_cfg$classification$primary_p_threshold)) {
    nominal_p_cut <- run_cfg$classification$primary_p_threshold
  }

  out <- lapply(genes, function(g) {
    a_all <- attrition_df %>% filter(gene == g)
    a_primary <- a_all %>% filter(outcome_name == "RA_EUR_PRIMARY")
    p1 <- preferred_df %>% filter(gene == g, outcome_name == "RA_EUR_PRIMARY")
    rp <- preferred_df %>% filter(gene == g, role == "replication") %>% arrange(outcome_name)
    sp <- preferred_df %>% filter(gene == g, outcome_name == "RA_EUR_SEROPOS")
    sn <- preferred_df %>% filter(gene == g, outcome_name == "RA_EUR_SERONEG")
    eas <- preferred_df %>% filter(gene == g, outcome_name == "RA_EAS_BBJ")
    tm <- if (!is.null(target_meta)) target_meta %>% filter(gene == g) else NULL

    primary_non_est <- nrow(a_primary) > 0 && all(!a_primary$estimable)
    primary_non_est_reason <- if (nrow(a_primary) > 0) {
      a_primary$non_estimable_reason[which(!is.na(a_primary$non_estimable_reason))[1]] %||% NA_character_
    } else {
      NA_character_
    }

    if (primary_non_est || nrow(p1) == 0) {
      return(data.frame(
        gene = g,
        class = "Non_estimable",
        support_subtype = NA_character_,
        primary_nominal_support = NA,
        replication_direction_concordant = NA,
        replication_nominal_support = NA,
        replication_concordant_n = NA_integer_,
        serostatus_direction_concordant = NA,
        eas_direction_concordant = NA,
        primary_p = NA_real_,
        replication_p = NA_real_,
        primary_beta = NA_real_,
        replication_beta = NA_real_,
        expected_direction_concordant = NA,
        alignment_caution_flag = NA_character_,
        non_estimable_reason = primary_non_est_reason,
        stringsAsFactors = FALSE
      ))
    }

    dir_mod <- expected_direction_modifier(tm, p1$b[1])
    primary_sig <- safe_sig(p1, nominal_p_cut)
    primary_supportive <- isTRUE(primary_sig) && !identical(dir_mod$expected_direction_concordant, FALSE)

    rep_same_sign_vec <- if (nrow(rp) > 0) vapply(rp$b, function(x) same_dir(p1$b[1], x), logical(1)) else logical(0)
    rep_sig_vec <- if (nrow(rp) > 0) vapply(seq_len(nrow(rp)), function(i) safe_sig(rp[i, , drop = FALSE], nominal_p_cut), logical(1)) else logical(0)
    rep_strong_vec <- if (length(rep_same_sign_vec) > 0 && length(rep_sig_vec) > 0) rep_same_sign_vec & rep_sig_vec else logical(0)

    rep_same_sign <- if (length(rep_same_sign_vec) > 0) any(rep_same_sign_vec, na.rm = TRUE) else NA
    rep_same_sign_n <- if (length(rep_same_sign_vec) > 0) sum(rep_same_sign_vec, na.rm = TRUE) else NA_integer_
    rep_nominal_sig <- if (length(rep_sig_vec) > 0) any(rep_sig_vec, na.rm = TRUE) else NA
    rep_strong_support <- if (length(rep_strong_vec) > 0) any(rep_strong_vec, na.rm = TRUE) else NA

    sero_same_sign <- nrow(sp) > 0 && nrow(sn) > 0 && same_dir(sp$b[1], sn$b[1])
    eas_same_sign <- nrow(eas) > 0 && same_dir(p1$b[1], eas$b[1])

    cls <- "Not_supported"
    subtype <- "null"

    if (isTRUE(primary_supportive) && isTRUE(rep_strong_support)) {
      cls <- "Strong"
      subtype <- NA_character_
    } else if (isTRUE(primary_supportive) && (isTRUE(rep_same_sign) || isTRUE(sero_same_sign) || isTRUE(eas_same_sign))) {
      cls <- "Weak_or_partial"
      subtype <- NA_character_
    } else if (isTRUE(primary_supportive)) {
      cls <- "Weak_or_partial"
      subtype <- NA_character_
    } else if (isTRUE(primary_sig) && identical(dir_mod$expected_direction_concordant, FALSE)) {
      cls <- "Not_supported"
      subtype <- "inverse"
    }

    data.frame(
      gene = g,
      class = cls,
      support_subtype = subtype,
      primary_nominal_support = primary_supportive,
      replication_direction_concordant = rep_same_sign,
      replication_nominal_support = rep_nominal_sig,
      replication_concordant_n = rep_same_sign_n,
      serostatus_direction_concordant = sero_same_sign,
      eas_direction_concordant = eas_same_sign,
      primary_p = p1$pval[1],
      replication_p = first_or_na(rp, "pval"),
      primary_beta = p1$b[1],
      replication_beta = first_or_na(rp, "b"),
      expected_direction_concordant = dir_mod$expected_direction_concordant,
      alignment_caution_flag = dir_mod$alignment_caution_flag,
      non_estimable_reason = NA_character_,
      stringsAsFactors = FALSE
    )
  })

  bind_rows(out)
}
