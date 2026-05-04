# R/16_failure_mode_cat6.R
# Category 6: post-transcriptional discordance (extended implementation)
# Version: 1.0
# Author: Mihye Kwon
#
# Dependencies:
#   - R/15_failure_mode_diagnostic.R
#   - configs/failure_mode_diagnostic.yaml
suppressPackageStartupMessages({
  library(dplyr)
  library(httr)
  library(jsonlite)
})

# ─────────────────────────────────────────────────────────────────────────────
# R/16_failure_mode_cat6.R
#
#
#   Rule 2  ARIES mQTL  — MRInstruments::aries_mqtl (aries_mqtl_17genes.csv)
#
#
# ─────────────────────────────────────────────────────────────────────────────


# ─────────────────────────────────────────────────────────────────────────────
# CONSTANTS
# ─────────────────────────────────────────────────────────────────────────────

# sQTL datasets (eQTL Catalogue API v2)
CAT6_SQTL_DATASETS <- c(
  "QTD000359",   # GTEx v8 whole blood txrev (primary — confirmed HTTP 200)
  "QTD000360",   # GTEx v8 whole blood leafcutter
  "QTD000025",   # BLUEPRINT monocyte leafcutter
  "QTD000035",   # BLUEPRINT CD4T leafcutter
  "QTD000478",   # Schmiedel 2018 B cell leafcutter
  "QTD000553"    # additional immune sQTL
)

CAT6_SQTL_KNOWN <- list(
  CAT6_SQTL_DISCORDANCE = c("cd86", "fcgr1a", "fcgr2b", "fcgr3b", "jak1", "jak2"),
  SPLICING_QTL_DETECTED  = c("fcgr2a", "jak3")
)

CAT6_PT_HIGH     <- c("cd80", "il1r1", "ms4a1", "tnf")
CAT6_PT_MODERATE <- c("lta")

ARIES_PVAL_THRESHOLD <- 1e-5

CAT6_RECOMMENDED_STRATEGY <- list(
  cd80   = "Multi-protein co-stimulation interaction MR or pQTL-MR (CD80/CD86/CTLA4)",
  cd86   = "Multi-protein co-stimulation interaction MR or pQTL-MR (CD80/CD86/CTLA4)",
  fcgr1a = "Pathway-level multi-target MR or restrict to direct drug target",
  fcgr2a = "Pathway-level multi-target MR or restrict to direct drug target",
  fcgr2b = "Pathway-level multi-target MR or restrict to direct drug target",
  fcgr2c = "Pathway-level multi-target MR or restrict to direct drug target",
  fcgr3a = "Pathway-level multi-target MR or restrict to direct drug target",
  fcgr3b = "Pathway-level multi-target MR or restrict to direct drug target",
  il1r1  = "pQTL-MR distinguishing membrane-bound vs soluble receptor",
  jak1   = "Coding variant MR (kinase domain missense) or pQTL-MR",
  jak2   = "Coding variant MR (kinase domain missense) or pQTL-MR",
  jak3   = "Coding variant MR (kinase domain missense) or pQTL-MR",
  lta    = "pQTL-MR using circulating protein (UKB-PPP, deCODE)",
  ms4a1  = "Immune cell abundance MR (B-cell GWAS) or scQTL-MR",
  tnf    = "pQTL-MR using circulating protein (UKB-PPP, deCODE)"
)

# Note: eQTL Catalogue v2 accepts unversioned IDs; versioned IDs used as per user spec
CAT6_ENSEMBL_IDS <- list(
  cd80   = "ENSG00000121594",
  cd86   = "ENSG00000114013",
  fcgr1a = "ENSG00000150337",
  fcgr2a = "ENSG00000143226",
  fcgr2b = "ENSG00000072694",
  fcgr2c = "ENSG00000162747",
  fcgr3a = "ENSG00000203747",
  fcgr3b = "ENSG00000162757",
  il1r1  = "ENSG00000115594",
  jak1   = "ENSG00000162434",
  jak2   = "ENSG00000096968",
  jak3   = "ENSG00000105639",
  lta    = "ENSG00000204490",
  ms4a1  = "ENSG00000156738",
  tnf    = "ENSG00000232810"
)


# ─────────────────────────────────────────────────────────────────────────────
# UTILITIES
# ─────────────────────────────────────────────────────────────────────────────

cat6_log <- function(msg, ...) {
  cat(sprintf("[cat6 %s] %s\n",
    format(Sys.time(), "%H:%M:%S"),
    sprintf(msg, ...)))
}

safe_csv <- function(path) {
  if (!file.exists(path)) return(data.frame())
  tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) { warning("[cat6] Cannot read ", path, ": ", conditionMessage(e)); data.frame() }
  )
}

write_csv_safe <- function(df, path) {
  tryCatch({
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(df, file = path, row.names = FALSE)
    cat6_log("Written: %s (%d rows)", basename(path), nrow(df))
  }, error = function(e) warning("[cat6] write_csv_safe: ", conditionMessage(e)))
  invisible(path)
}


# ─────────────────────────────────────────────────────────────────────────────
# RULE 1: sQTL via eQTL Catalogue API v2
#
#
# ─────────────────────────────────────────────────────────────────────────────

query_sqtl_eqtl_catalogue <- function(target_genes, cache_dir, delay_ms = 500) {

  eqtl_base <- "https://www.ebi.ac.uk/eqtl/api/v2"

  query_one <- function(gene, dataset_id) {
    ensembl_id <- CAT6_ENSEMBL_IDS[[gene]]
    if (is.null(ensembl_id)) return(NULL)

    cache_key  <- file.path(cache_dir, sprintf("sqtl_%s_%s.rds", gene, dataset_id))
    if (file.exists(cache_key)) {
      cached <- tryCatch(readRDS(cache_key), error = function(e) NULL)
      if (!is.null(cached)) return(cached)
    }

    Sys.sleep(delay_ms / 1000)
    url <- sprintf("%s/datasets/%s/associations?gene_id=%s&nlog10p=5&size=10",
                   eqtl_base, dataset_id, ensembl_id)
    resp <- tryCatch(
      httr::GET(url, httr::timeout(30), httr::accept_json()),
      error = function(e) NULL
    )

    result <- NULL
    if (!is.null(resp) && httr::status_code(resp) == 200L) {
      parsed <- tryCatch(
        httr::content(resp, as = "parsed", type = "application/json", encoding = "UTF-8"),
        error = function(e) NULL
      )
      # v2 API: raw array (not resp$results)
      if (is.null(parsed)) {
        result <- list()
      } else if (is.list(parsed) && is.null(names(parsed))) {
        result <- parsed       # plain array → hits
      } else {
        result <- parsed$results %||% list()
      }
    }
    # HTTP 400 = No results; NULL = network error → return NULL for retry detection
    if (!is.null(result)) {
      tryCatch({
        dir.create(dirname(cache_key), showWarnings = FALSE, recursive = TRUE)
        saveRDS(result, cache_key)
      }, error = function(e) invisible(NULL))
    }
    result
  }

  result <- lapply(target_genes, function(g) {
    cat6_log("sQTL query: %s", g)
    api_worked   <- FALSE
    hit_datasets <- character(0)

    for (ds in CAT6_SQTL_DATASETS) {
      res <- query_one(g, ds)
      if (!is.null(res)) {
        api_worked <- TRUE
        if (length(res) > 0) hit_datasets <- c(hit_datasets, ds)
      }
    }

    has_sqtl <- length(hit_datasets) > 0

    sqtl_source <- if (api_worked) "eqtl_catalogue_api_v2" else "known_results_fallback"
    if (!api_worked) {
      all_known <- c(CAT6_SQTL_KNOWN$CAT6_SQTL_DISCORDANCE, CAT6_SQTL_KNOWN$SPLICING_QTL_DETECTED)
      has_sqtl  <- g %in% all_known
      if (has_sqtl) hit_datasets <- "fallback"
    }

    list(
      gene          = g,
      has_sqtl      = has_sqtl,
      sqtl_datasets = paste(hit_datasets, collapse = ";"),
      sqtl_source   = sqtl_source,
      api_worked    = api_worked
    )
  })
  names(result) <- target_genes
  result
}


# ─────────────────────────────────────────────────────────────────────────────
# RULE 2: ARIES mQTL
#
# ─────────────────────────────────────────────────────────────────────────────

load_aries_mqtl <- function(results_dir, cache_dir) {
  csv_path <- file.path(results_dir, "aries_mqtl_17genes.csv")
  if (!file.exists(csv_path)) {
    cat6_log("aries_mqtl_17genes.csv not found in results_dir, checking cache/")
    csv_path <- NULL
  }

  if (!is.null(csv_path)) {
    aries_df <- safe_csv(csv_path)
    cat6_log("ARIES CSV loaded: %d rows from %s", nrow(aries_df), basename(csv_path))
    return(summarise_aries_csv(aries_df))
  }

  cat6_log("Falling back to individual aries_mqtl_*.rds cache files")
  rds_files <- list.files(cache_dir, pattern = "^aries_mqtl_.+\\.rds$", full.names = TRUE)
  if (length(rds_files) == 0) {
    cat6_log("No ARIES mQTL cache found")
    return(list())
  }

  all_rows <- lapply(rds_files, function(f) {
    df <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)
    df
  })
  combined <- do.call(rbind, Filter(Negate(is.null), all_rows))
  if (is.null(combined) || nrow(combined) == 0) return(list())
  summarise_aries_csv(combined)
}


# ARIES CSV/RDS → per-gene summary
summarise_aries_csv <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(list())

  if ("gene" %in% names(df)) {
    df$gene_lc <- tolower(df$gene)
  } else {
    return(list())
  }

  genes <- unique(df$gene_lc)
  result <- lapply(genes, function(g) {
    sub <- df[df$gene_lc == g, , drop = FALSE]

    if ("pval" %in% names(sub)) {
      pvs  <- suppressWarnings(as.numeric(sub$pval))
      best_idx <- which.min(pvs)
    } else {
      best_idx <- 1L
      pvs <- rep(NA_real_, nrow(sub))
    }

    best_pval <- pvs[best_idx]

    tested_positive <- any(!is.na(pvs) & pvs < ARIES_PVAL_THRESHOLD)
    mqtl_status <- if (tested_positive) "tested_positive" else "tested_negative"

    # best CpG
    best_cpg <- if ("cpg" %in% names(sub)) as.character(sub$cpg[best_idx]) else NA_character_

    # gene_location at best CpG
    gene_loc <- if ("gene_location" %in% names(sub)) as.character(sub$gene_location[best_idx]) else NA_character_

    sig_sub <- sub[!is.na(pvs) & pvs < ARIES_PVAL_THRESHOLD, , drop = FALSE]
    n_cpg       <- if ("cpg" %in% names(sig_sub) && nrow(sig_sub) > 0)
      length(unique(sig_sub$cpg)) else NA_integer_
    n_timepoints <- if ("timepoint" %in% names(sig_sub) && nrow(sig_sub) > 0)
      length(unique(sig_sub$timepoint)) else NA_integer_

    list(
      gene          = g,
      mqtl_status   = mqtl_status,
      aries_best_pval = best_pval,
      aries_best_cpg  = best_cpg,
      aries_n_cpg     = n_cpg,
      aries_n_timepoints = n_timepoints,
      aries_gene_location = gene_loc,
      source         = "aries_mrinstruments"
    )
  })
  names(result) <- genes
  result
}


# ─────────────────────────────────────────────────────────────────────────────
#
# cd80, il1r1, ms4a1, tnf → High
# ─────────────────────────────────────────────────────────────────────────────

get_pt_mismatch_flag <- function(gene) {
  g <- tolower(gene)
  if (g %in% CAT6_PT_HIGH)     return(list(flag = "POST_TRANSCRIPTIONAL_MISMATCH", confidence = "High"))
  if (g %in% CAT6_PT_MODERATE) return(list(flag = "POST_TRANSCRIPTIONAL_MISMATCH", confidence = "Moderate"))
  list(flag = NA_character_, confidence = NA_character_)
}


# ─────────────────────────────────────────────────────────────────────────────
#
# ─────────────────────────────────────────────────────────────────────────────

validate_eqtl_pqtl_discordant <- function(gene, current_flag, pqtl_support_df) {
  if (is.null(pqtl_support_df) || nrow(pqtl_support_df) == 0) return(current_flag)
  if (!identical(current_flag, "EQTL_PQTL_DISCORDANT")) return(current_flag)

  g <- tolower(gene)
  row <- pqtl_support_df[tolower(pqtl_support_df$gene) == g, , drop = FALSE]
  if (nrow(row) == 0) return(current_flag)

  n_primary <- NA_integer_
  if ("pqtl_n_sources_primary_support" %in% names(row)) {
    n_primary <- suppressWarnings(as.integer(row$pqtl_n_sources_primary_support[1]))
  } else if (ncol(row) >= 7) {
    n_primary <- suppressWarnings(as.integer(row[[7]][1]))
  }

  if (!is.na(n_primary) && n_primary == 0L) {
    cat6_log("Rule 4: %s EQTL_PQTL_DISCORDANT → POST_TRANSCRIPTIONAL_MISMATCH (pqtl_n_sources_primary_support=0)", g)
    return("POST_TRANSCRIPTIONAL_MISMATCH_corrected")
  }
  current_flag
}


# ─────────────────────────────────────────────────────────────────────────────
# CAT6 FLAG PRIORITY
#
# Priority: EQTL_PQTL_DISCORDANT > CAT6_SQTL_DISCORDANCE > SPLICING_QTL_DETECTED
#           > ALTERNATIVE_QTL_AVAILABLE > POST_TRANSCRIPTIONAL_MISMATCH
#
# sQTL vs eQTL estimability:
# ─────────────────────────────────────────────────────────────────────────────

assign_cat6_flag <- function(gene, sqtl_info, aries_info, pt_flag, blood_eqtl_estimable) {
  g <- tolower(gene)

  has_sqtl  <- isTRUE(sqtl_info$has_sqtl)
  has_mqtl  <- !is.null(aries_info) && identical(aries_info$mqtl_status, "tested_positive")
  has_blood_eqtl <- isTRUE(blood_eqtl_estimable)

  if (has_sqtl && has_blood_eqtl)  return(list(flag = "CAT6_SQTL_DISCORDANCE",      confidence = "Moderate"))
  if (has_sqtl && !has_blood_eqtl) return(list(flag = "SPLICING_QTL_DETECTED",       confidence = "Moderate"))
  if (has_mqtl && !has_blood_eqtl) return(list(flag = "ALTERNATIVE_QTL_AVAILABLE",   confidence = "Moderate"))
  if (!is.na(pt_flag$flag))         return(list(flag = pt_flag$flag, confidence = pt_flag$confidence))

  list(flag = NA_character_, confidence = NA_character_)
}


# ─────────────────────────────────────────────────────────────────────────────
# MAIN: run_cat6_full
#
#   config_path  — failure_mode_diagnostic.yaml
# ─────────────────────────────────────────────────────────────────────────────

run_cat6_full <- function(results_dir,
                          config_path  = "configs/failure_mode_diagnostic.yaml",
                          cache_dir    = "cache/failmode_api",
                          output_dir   = results_dir,
                          run_sqtl_api = FALSE) {

  cat6_log("START results_dir=%s run_sqtl_api=%s", results_dir, run_sqtl_api)

  comp_path <- file.path(results_dir, "fail_diagnostic_composite.csv")
  composite <- safe_csv(comp_path)
  if (nrow(composite) == 0) {
    cat6_log("ERROR: fail_diagnostic_composite.csv not found or empty at %s", comp_path)
    return(invisible(NULL))
  }
  cat6_log("Composite loaded: %d genes x %d cols", nrow(composite), ncol(composite))

  target_genes <- tolower(composite$gene)

  pqtl_support_df <- safe_csv(file.path(results_dir, "pqtl_support_summary.csv"))
  if (nrow(pqtl_support_df) > 0) {
    if (!"gene" %in% names(pqtl_support_df)) {
      names(pqtl_support_df)[1] <- "gene"
      if (ncol(pqtl_support_df) >= 7) names(pqtl_support_df)[7] <- "pqtl_n_sources_primary_support"
    }
    cat6_log("pqtl_support_summary loaded: %d rows", nrow(pqtl_support_df))
  }

  # ── 3. Rule 1: sQTL ───────────────────────────────────────────────────────
  if (isTRUE(run_sqtl_api)) {
    cat6_log("Rule 1: querying eQTL Catalogue API v2 for sQTL")
    sqtl_results <- query_sqtl_eqtl_catalogue(target_genes, cache_dir)
  } else {
    cat6_log("Rule 1: using hard-coded known sQTL results (run_sqtl_api=FALSE)")
    sqtl_results <- lapply(target_genes, function(g) {
      if (g %in% CAT6_SQTL_KNOWN$CAT6_SQTL_DISCORDANCE) {
        list(gene = g, has_sqtl = TRUE, sqtl_datasets = "known", sqtl_source = "known_results", api_worked = FALSE)
      } else if (g %in% CAT6_SQTL_KNOWN$SPLICING_QTL_DETECTED) {
        list(gene = g, has_sqtl = TRUE, sqtl_datasets = "known", sqtl_source = "known_results", api_worked = FALSE)
      } else {
        list(gene = g, has_sqtl = FALSE, sqtl_datasets = NA_character_, sqtl_source = "known_results", api_worked = FALSE)
      }
    })
    names(sqtl_results) <- target_genes
  }

  # ── 4. Rule 2: ARIES mQTL ─────────────────────────────────────────────────
  cat6_log("Rule 2: loading ARIES mQTL data")
  aries_results <- load_aries_mqtl(results_dir, cache_dir)
  # genes not found in ARIES → tested_negative
  for (g in target_genes) {
    if (is.null(aries_results[[g]])) {
      aries_results[[g]] <- list(
        gene = g, mqtl_status = "tested_negative",
        aries_best_pval = NA_real_, aries_best_cpg = NA_character_,
        aries_n_cpg = NA_integer_, aries_n_timepoints = NA_integer_,
        aries_gene_location = NA_character_, source = "aries_mrinstruments"
      )
    }
  }

  # ── 5. blood eQTL estimability (attrition cascade: S5 == TRUE) ────────────
  blood_eqtl_map <- setNames(
    composite$class != "Non_estimable",
    tolower(composite$gene)
  )

  detail_rows <- lapply(target_genes, function(g) {
    sqtl_info <- sqtl_results[[g]] %||% list(has_sqtl = FALSE)
    aries_info <- aries_results[[g]]
    pt_flag    <- get_pt_mismatch_flag(g)
    blood_eqtl <- isTRUE(blood_eqtl_map[[g]])

    flag_result <- assign_cat6_flag(g, sqtl_info, aries_info, pt_flag, blood_eqtl)

    current_flag <- composite$cat6_flag[tolower(composite$gene) == g]
    current_flag <- if (length(current_flag) > 0 && !is.na(current_flag[1])) current_flag[1] else NA_character_

    validated_flag <- validate_eqtl_pqtl_discordant(g, current_flag, pqtl_support_df)
    if (identical(validated_flag, "POST_TRANSCRIPTIONAL_MISMATCH_corrected")) {
      flag_result <- list(flag = "POST_TRANSCRIPTIONAL_MISMATCH", confidence = "Moderate")
    }

    mqtl_pos <- identical(aries_info$mqtl_status, "tested_positive")
    aries_flag       <- if (mqtl_pos) "ALTERNATIVE_QTL_AVAILABLE" else NA_character_
    aries_confidence <- if (mqtl_pos) "Moderate" else NA_character_
    aries_source     <- if (mqtl_pos) "aries_mrinstruments" else NA_character_

    data.frame(
      gene                    = g,
      cat6_flag               = flag_result$flag,
      cat6_confidence         = flag_result$confidence,
      cat6_aries_flag         = aries_flag,
      cat6_aries_confidence   = aries_confidence,
      cat6_aries_source       = aries_source,
      mqtl_status             = aries_info$mqtl_status %||% "data_unavailable",
      aries_best_pval         = aries_info$aries_best_pval %||% NA_real_,
      aries_n_cpg             = aries_info$aries_n_cpg    %||% NA_integer_,
      aries_n_timepoints      = aries_info$aries_n_timepoints %||% NA_integer_,
      aries_gene_location     = aries_info$aries_gene_location %||% NA_character_,
      blood_eqtl_estimable    = blood_eqtl,
      sqtl_detected           = isTRUE(sqtl_info$has_sqtl),
      sqtl_datasets           = sqtl_info$sqtl_datasets %||% NA_character_,
      sqtl_source             = sqtl_info$sqtl_source   %||% NA_character_,
      rule4_corrected         = identical(validated_flag, "POST_TRANSCRIPTIONAL_MISMATCH_corrected"),
      annotation_timestamp    = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
      stringsAsFactors = FALSE
    )
  })
  cat6_detail <- dplyr::bind_rows(detail_rows)
  n_flagged <- sum(!is.na(cat6_detail$cat6_flag))
  cat6_log("CAT6 flags assigned: %d/%d genes flagged", n_flagged, nrow(cat6_detail))

  updated <- update_composite(composite, cat6_detail)
  cat6_log("Composite updated: %d rows x %d cols", nrow(updated), ncol(updated))

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  write_csv_safe(cat6_detail, file.path(output_dir, "cat6_full_detail.csv"))
  write_csv_safe(updated,     file.path(output_dir, "fail_diagnostic_composite.csv"))

  cat6_log("DONE")
  invisible(list(detail = cat6_detail, composite = updated))
}


# ─────────────────────────────────────────────────────────────────────────────
# UPDATE_COMPOSITE
#
# ─────────────────────────────────────────────────────────────────────────────

update_composite <- function(composite, cat6_detail) {

  composite$gene_lc <- tolower(composite$gene)

  for (i in seq_len(nrow(composite))) {
    g   <- composite$gene_lc[i]
    row <- cat6_detail[tolower(cat6_detail$gene) == g, , drop = FALSE]
    if (nrow(row) == 0) next
    composite$cat6_flag[i]       <- row$cat6_flag[1]
    composite$cat6_confidence[i] <- row$cat6_confidence[1]
  }

  aries_cols <- c("cat6_aries_flag", "cat6_aries_confidence", "cat6_aries_source",
                  "mqtl_status", "aries_best_pval", "aries_n_cpg",
                  "aries_n_timepoints", "aries_gene_location")
  for (col in aries_cols) {
    composite[[col]] <- NA
  }
  for (i in seq_len(nrow(composite))) {
    g   <- composite$gene_lc[i]
    row <- cat6_detail[tolower(cat6_detail$gene) == g, , drop = FALSE]
    if (nrow(row) == 0) next
    for (col in aries_cols) {
      if (col %in% names(row)) composite[[col]][i] <- row[[col]][1]
    }
  }

  for (i in seq_len(nrow(composite))) {
    g        <- composite$gene_lc[i]
    strategy <- CAT6_RECOMMENDED_STRATEGY[[g]]
    if (!is.null(strategy)) composite$recommended_strategy[i] <- strategy
  }

  composite$n_categories_flagged <- apply(composite, 1, function(row) {
    flags_per_cat <- c(
      cat1 = !is.na(row["cat1_flag"]),
      cat2 = !is.na(row["cat2_flag"]),
      cat3 = !is.na(row["cat3_flag"]),
      cat4 = !is.na(row["cat4_flag"]),
      cat5 = !is.na(row["cat5_flag"]),
      cat6 = (!is.na(row["cat6_flag"])) | (!is.na(row["cat6_aries_flag"])),
      cat7 = !is.na(row["cat7_flag"])
    )
    sum(as.logical(flags_per_cat), na.rm = TRUE)
  })

  ordered_cols <- c(
    "gene", "class",
    "cat1_flag", "cat1_confidence",
    "cat2_flag", "cat2_confidence", "cat2_locus_score",
    "cat3_flag", "cat3_confidence",
    "cat4_flag", "cat4_confidence",
    "cat5_flag", "cat5_confidence",
    "cat6_flag", "cat6_confidence",
    "cat6_aries_flag", "cat6_aries_confidence", "cat6_aries_source",
    "mqtl_status", "aries_best_pval", "aries_n_cpg",
    "aries_n_timepoints", "aries_gene_location",
    "cat7_flag", "cat7_confidence",
    "n_categories_flagged", "recommended_strategy", "annotation_timestamp"
  )

  for (col in ordered_cols) {
    if (!col %in% names(composite)) composite[[col]] <- NA
  }

  extra_cols <- setdiff(names(composite), c(ordered_cols, "gene_lc"))
  final_cols <- c(ordered_cols, extra_cols)
  composite <- composite[, final_cols, drop = FALSE]

  composite$annotation_timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  composite
}


# ─────────────────────────────────────────────────────────────────────────────
# EXECUTION BLOCK
# ─────────────────────────────────────────────────────────────────────────────

if (!exists(".cat6_sourced_only", mode = "logical")) {

  run_dir <- tryCatch({
    if (exists("results_dir")) {
      results_dir
    } else {
      run_parent <- "runs"
      run_dirs   <- list.dirs(run_parent, recursive = FALSE, full.names = TRUE)
      if (length(run_dirs) == 0) stop("runs/ directory empty")
      latest_run <- tail(sort(run_dirs), 1)
      file.path(latest_run, "results")
    }
  }, error = function(e) {
    stop("[cat6] Cannot determine results_dir: ", conditionMessage(e))
  })

  cache_dir_use <- "cache/failmode_api"
  config_path   <- "configs/failure_mode_diagnostic.yaml"

  cat6_log("Executing run_cat6_full()")
  cat6_log("  results_dir : %s", run_dir)
  cat6_log("  cache_dir   : %s", cache_dir_use)
  cat6_log("  run_sqtl_api: FALSE (using known results)")

  cat6_out <- run_cat6_full(
    results_dir  = run_dir,
    config_path  = config_path,
    cache_dir    = cache_dir_use,
    output_dir   = run_dir,
  )

  if (!is.null(cat6_out)) {
    cat6_log("Summary:")
    cat6_log("  Composite rows  : %d", nrow(cat6_out$composite))
    cat6_log("  CAT6 flags set  : %d", sum(!is.na(cat6_out$detail$cat6_flag)))
    cat6_log("  ARIES mQTL pos  : %d", sum(cat6_out$detail$mqtl_status == "tested_positive", na.rm=TRUE))
    cat6_log("  Rule4 corrected : %d", sum(cat6_out$detail$rule4_corrected, na.rm=TRUE))
  }
}
