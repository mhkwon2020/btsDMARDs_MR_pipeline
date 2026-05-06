# R/03_instrument_builder.R
# Instrument builder utilities
# Version: 1.0
# Author: Mihye Kwon
suppressPackageStartupMessages({
  library(dplyr)
  library(TwoSampleMR)
  library(fs)
})

filter_target_variants <- function(exp_raw, target_gene) {
  exp_raw %>%
    dplyr::filter(.data$gene == tolower(.env$target_gene))
}

apply_p_threshold <- function(gdat, p_threshold) {
  gdat %>% filter(!is.na(pval) & pval <= p_threshold)
}

run_ld_clumping <- function(exp_dat, protocol_cfg, ancestry) {

  ref_bfile <- ifelse(
    toupper(ancestry) == "EAS",
    Sys.getenv("LD_REF_EAS"),
    Sys.getenv("LD_REF_EUR")
  )

  TwoSampleMR::clump_data(
    exp_dat,
    clump_r2 = protocol_cfg$instrument_selection$clumping$r2,
    clump_kb = protocol_cfg$instrument_selection$clumping$kb,
    clump_p1 = protocol_cfg$instrument_selection$p_threshold$primary,
    bfile = ref_bfile,
    plink_bin = Sys.getenv("PLINK_BIN")
  )
}

is_retryable_clump_error <- function(msg) {
  if (is.null(msg) || length(msg) == 0 || is.na(msg)) return(FALSE)
  msg_l <- tolower(paste(msg, collapse = " "))
  patterns <- c(
    "timed out",
    "timeout",
    "server code: 502",
    "server code: 503",
    "server code: 504",
    "bad gateway",
    "gateway",
    "temporarily blocked",
    "allowance",
    "connection",
    "could not resolve host",
    "failed to connect",
    "proxy",
    "opengwas exceeded 300 seconds"
  )
  any(vapply(patterns, grepl, logical(1), x = msg_l, fixed = TRUE))
}

get_clump_retry_attempts <- function(protocol_cfg) {
  val <- tryCatch(protocol_cfg$instrument_selection$clumping$retry$max_attempts, error = function(e) NULL)
  if (is.null(val) || !is.numeric(val) || length(val) != 1 || is.na(val) || val < 1) return(5L)
  as.integer(val)
}

get_clump_retry_delays <- function(protocol_cfg) {
  val <- tryCatch(protocol_cfg$instrument_selection$clumping$retry$delays_seconds, error = function(e) NULL)
  if (is.null(val) || !is.numeric(val) || length(val) == 0 || any(is.na(val))) {
    return(c(30, 60, 120, 300))
  }
  as.numeric(val)
}

get_clump_cache_dir <- function(protocol_cfg) {
  cache_dir <- tryCatch(protocol_cfg$environment$cache_dir, error = function(e) NULL)
  if (is.null(cache_dir) || !nzchar(cache_dir)) {
    cache_dir <- "cache/clumped_instruments"
  }
  cache_dir
}

make_clump_cache_path <- function(protocol_cfg, gene, ancestry) {
  cache_dir <- get_clump_cache_dir(protocol_cfg)
  fs::dir_create(cache_dir, recurse = TRUE)
  pth <- protocol_cfg$instrument_selection$p_threshold$primary
  r2 <- protocol_cfg$instrument_selection$clumping$r2
  kb <- protocol_cfg$instrument_selection$clumping$kb

  safe_gene <- gsub("[^A-Za-z0-9_\\-]", "_", tolower(gene))
  safe_anc <- gsub("[^A-Za-z0-9_\\-]", "_", toupper(ancestry))
  safe_p <- gsub("[^A-Za-z0-9_\\-]", "_", format(pth, scientific = TRUE))
  safe_r2 <- gsub("[^A-Za-z0-9_\\-]", "_", as.character(r2))
  safe_kb <- gsub("[^A-Za-z0-9_\\-]", "_", as.character(kb))

  fs::path(cache_dir, paste0(safe_gene, "__", safe_anc, "__p", safe_p, "__r2_", safe_r2, "__kb_", safe_kb, ".rds"))
}

run_ld_clumping_with_retry <- function(exp_dat, protocol_cfg, ancestry, gene,
                                       log_file = NULL, err_file = NULL) {
  max_attempts <- get_clump_retry_attempts(protocol_cfg)
  delays <- get_clump_retry_delays(protocol_cfg)

  last_error <- NULL

  for (attempt in seq_len(max_attempts)) {
    if (!is.null(log_file)) {
      log_message(log_file, "CLUMP TRY | gene={gene} ancestry={ancestry} attempt={attempt}/{max_attempts}")
    }

    res <- tryCatch(
      run_ld_clumping(exp_dat, protocol_cfg, ancestry),
      error = function(e) e
    )

    if (!inherits(res, "error")) {
      if (!is.null(log_file)) {
        n_clumped <- if (is.null(res)) 0L else nrow(res)
        log_message(log_file, "CLUMP OK | gene={gene} ancestry={ancestry} attempt={attempt}/{max_attempts} n={n_clumped}")
      }
      return(res)
    }

    last_error <- res
    msg <- conditionMessage(res)

    if (!is.null(err_file)) {
      log_message(err_file, "CLUMP FAIL | gene={gene} ancestry={ancestry} attempt={attempt}/{max_attempts} msg={msg}")
    }

    retryable <- is_retryable_clump_error(msg)

    if (!retryable) {
      stop(sprintf("Non-retryable clumping failure for gene=%s ancestry=%s: %s", gene, ancestry, msg), call. = FALSE)
    }

    if (attempt < max_attempts) {
      delay <- delays[min(attempt, length(delays))]
      if (!is.null(log_file)) {
        log_message(log_file, "CLUMP RETRY | gene={gene} ancestry={ancestry} sleeping_seconds={delay}")
      }
      Sys.sleep(delay)
    }
  }

  stop(
    sprintf(
      "Technical clumping failure after %d attempts for gene=%s ancestry=%s: %s",
      max_attempts, gene, ancestry, conditionMessage(last_error)
    ),
    call. = FALSE
  )
}

build_instrument_for_gene <- function(exp_raw, gene, ancestry, protocol_cfg, target_meta = NULL,
                                      log_file = NULL, err_file = NULL) {
  pth <- protocol_cfg$instrument_selection$p_threshold$primary

  gdat_raw <- filter_target_variants(exp_raw, gene)
  raw_n <- nrow(gdat_raw)

  gene_n_median <- if (raw_n > 0) unique(gdat_raw$n_eqtl_gene_median)[1] else NA_real_
  gene_n_max <- if (raw_n > 0) unique(gdat_raw$n_eqtl_gene_max)[1] else NA_real_
  gene_n_min <- if (raw_n > 0) unique(gdat_raw$n_eqtl_gene_min)[1] else NA_real_
  gene_n_distinct <- if (raw_n > 0) unique(gdat_raw$n_eqtl_gene_distinct)[1] else NA_real_

  base_row <- data.frame(
    gene = tolower(gene),
    ancestry = ancestry,
    n_eqtl_gene_median = gene_n_median,
    n_eqtl_gene_max = gene_n_max,
    n_eqtl_gene_min = gene_n_min,
    n_eqtl_gene_distinct = gene_n_distinct,
    raw_variant_count = raw_n,
    post_p_threshold_count = NA_integer_,
    post_clump_count = NA_integer_,
    post_outcome_extraction_count = NA_integer_,
    post_harmonisation_row_count = NA_integer_,
    post_harmonisation_usable_count = NA_integer_,
    palindromic_loss_count = NA_integer_,
    estimable = FALSE,
    non_estimable_reason = NA_character_,
    stringsAsFactors = FALSE
  )

  if (raw_n == 0) {
    base_row$non_estimable_reason <- "no_variant_after_target_filter"
    return(list(exp_dat = NULL, attrition = base_row))
  }

  gdat_thr <- apply_p_threshold(gdat_raw, pth)
  base_row$post_p_threshold_count <- nrow(gdat_thr)

  if (nrow(gdat_thr) == 0) {
    base_row$non_estimable_reason <- "no_variant_after_p_threshold"
    return(list(exp_dat = NULL, attrition = base_row))
  }

  gdat_thr <- gdat_thr %>%
    arrange(SNP, pval) %>%
    distinct(SNP, .keep_all = TRUE)

  if (!"n_eqtl_variant" %in% names(gdat_thr)) {
    if ("n_eqtl_gene_median" %in% names(gdat_thr)) {
      gdat_thr$n_eqtl_variant <- gdat_thr$n_eqtl_gene_median
    } else {
      gdat_thr$n_eqtl_variant <- NA_real_
    }
  }

  gdat_thr_df <- as.data.frame(gdat_thr)

  exp_dat <- TwoSampleMR::format_data(
    gdat_thr_df,
    type = "exposure",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    chr_col = "chr",
    pos_col = "pos"
  )

  exp_dat$exposure <- tolower(gene)
  exp_dat$id.exposure <- tolower(gene)

  meta_cols <- gdat_thr_df %>%
    dplyr::select(
      SNP,
      n_eqtl_variant,
      n_eqtl_gene_median,
      n_eqtl_gene_max
    ) %>%
    dplyr::distinct()

  exp_dat <- exp_dat %>%
    left_join(meta_cols, by = c("SNP" = "SNP"))

  cache_path <- make_clump_cache_path(protocol_cfg, gene, ancestry)

  clumped <- NULL
  if (file.exists(cache_path)) {
    if (!is.null(log_file)) {
      log_message(log_file, "CLUMP CACHE HIT | gene={gene} ancestry={ancestry} path={cache_path}")
    }
    clumped <- readRDS(cache_path)
  } else {
    if (!is.null(log_file)) {
      log_message(log_file, "CLUMP CACHE MISS | gene={gene} ancestry={ancestry} path={cache_path}")
    }
    clumped <- run_ld_clumping_with_retry(
      exp_dat = exp_dat,
      protocol_cfg = protocol_cfg,
      ancestry = ancestry,
      gene = gene,
      log_file = log_file,
      err_file = err_file
    )
    saveRDS(clumped, cache_path)
    if (!is.null(log_file)) {
      log_message(log_file, "CLUMP CACHE WRITE | gene={gene} ancestry={ancestry} path={cache_path}")
    }
  }

  base_row$post_clump_count <- if (is.null(clumped)) 0L else length(unique(clumped$SNP))

  if (is.null(clumped) || nrow(clumped) == 0) {
    base_row$non_estimable_reason <- "no_variant_after_clumping"
    return(list(exp_dat = NULL, attrition = base_row))
  }

  list(exp_dat = clumped, attrition = base_row)
}
