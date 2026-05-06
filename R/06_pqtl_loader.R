# R/06_pqtl_loader.R
# pQTL loader utilities
# Version: 1.0
# Author: Mihye Kwon
suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
})

validate_pqtl_manifest <- function(manifest, target_genes = NULL) {
  check_nonempty(manifest$version, "pqtl_manifest.version")
  check_nonempty(manifest$sources, "pqtl_manifest.sources")

  required_cols <- c(
    "gene", "snp", "beta", "se", "pval", "eaf",
    "chr", "pos", "effect_allele", "other_allele"
  )

  for (src in manifest$sources) {
    check_nonempty(src$source_name, "pqtl_manifest.sources[[i]]$source_name")
    check_true_false(src$enabled, paste0("pqtl_manifest source enabled flag for ", src$source_name))
    if (isTRUE(src$enabled)) {
      check_nonempty(src$file, paste0("pqtl_manifest file for ", src$source_name))
      check_nonempty(src$colmap, paste0("pqtl_manifest colmap for ", src$source_name))
      missing_map <- setdiff(required_cols, names(src$colmap))
      if (length(missing_map) > 0) {
        fail("pQTL source %s is missing required colmap entries: %s",
             src$source_name, paste(missing_map, collapse = ", "))
      }
    }
  }

  invisible(TRUE)
}

read_pqtl_manifest <- function(path, strict = TRUE) {
  if (!file.exists(path)) {
    if (isTRUE(strict)) {
      fail("pQTL source manifest does not exist: %s", path)
    }
    return(list(version = "missing_manifest", sources = list()))
  }

  manifest <- yaml::read_yaml(path)
  if (is.null(manifest)) {
    if (isTRUE(strict)) {
      fail("pQTL source manifest could not be parsed: %s", path)
    }
    return(list(version = "empty_manifest", sources = list()))
  }

  manifest
}

standardize_pqtl_source <- function(raw_df, source_cfg, target_genes = NULL) {
  cm <- source_cfg$colmap

  pull_col <- function(nm, required = TRUE, default = NA) {
    col_name <- cm[[nm]] %||% NULL
    if (is.null(col_name)) {
      if (required) fail("Required column mapping '%s' missing for %s", nm, source_cfg$source_name)
      return(rep(default, nrow(raw_df)))
    }
    if (!col_name %in% names(raw_df)) {
      if (required) fail("Mapped column '%s' not found in %s for source %s",
                         col_name, nm, source_cfg$source_name)
      return(rep(default, nrow(raw_df)))
    }
    raw_df[[col_name]]
  }

  out <- data.frame(
    gene = tolower(as.character(pull_col("gene"))),
    SNP = as.character(pull_col("snp")),
    beta = as.numeric(pull_col("beta")),
    se = as.numeric(pull_col("se")),
    pval = as.numeric(pull_col("pval")),
    eaf = as.numeric(pull_col("eaf")),
    chr = suppressWarnings(as.numeric(pull_col("chr"))),
    pos = suppressWarnings(as.numeric(pull_col("pos"))),
    effect_allele = as.character(pull_col("effect_allele")),
    other_allele = as.character(pull_col("other_allele")),
    stringsAsFactors = FALSE
  )

  out$maf <- pmin(out$eaf, 1 - out$eaf)
  out$pqtl_source <- source_cfg$source_name
  out$pqtl_source_id <- source_cfg$source_id %||% source_cfg$source_name
  out$proxy_layer <- "pQTL"
  out$n_eqtl_variant <- as.numeric(source_cfg$sample_size %||% NA_real_)
  out$n_eqtl_gene_median <- as.numeric(source_cfg$sample_size %||% NA_real_)
  out$n_eqtl_gene_max <- as.numeric(source_cfg$sample_size %||% NA_real_)
  out$n_eqtl_gene_min <- as.numeric(source_cfg$sample_size %||% NA_real_)
  out$n_eqtl_gene_distinct <- 1L

  out <- out %>%
    filter(
      !is.na(gene),
      nzchar(gene),
      !is.na(SNP),
      nzchar(SNP),
      !is.na(beta),
      !is.na(se),
      !is.na(pval),
      !is.na(eaf),
      !is.na(chr),
      !is.na(pos),
      !is.na(effect_allele),
      nzchar(effect_allele),
      !is.na(other_allele),
      nzchar(other_allele)
    ) %>%
    arrange(gene, pval) %>%
    distinct(gene, SNP, .keep_all = TRUE)

  if (!is.null(target_genes)) {
    out <- out %>% filter(gene %in% tolower(target_genes))
  }

  out
}

load_single_pqtl_source <- function(source_cfg, target_genes = NULL) {
  if (!isTRUE(source_cfg$enabled)) return(NULL)

  # Resolve ${VAR} placeholders in the file path (see 00_utils.R::expand_env_path)
  file_path <- expand_env_path(source_cfg$file)

  if (!file.exists(file_path)) {
    stop(
      "pQTL source file not found for '", source_cfg$source_name, "': ", file_path,
      "\nCheck that PQTL_DATA_DIR is set correctly in your .Renviron or environment.",
      call. = FALSE
    )
  }

  delim <- source_cfg$delimiter %||% "\t"
  raw_df <- utils::read.table(
    file_path,
    header = TRUE,
    sep = delim,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )

  standardize_pqtl_source(raw_df, source_cfg, target_genes = target_genes)
}

load_all_pqtl_sources <- function(manifest, target_genes = NULL) {
  enabled_sources <- Filter(function(x) isTRUE(x$enabled), manifest$sources %||% list())
  if (length(enabled_sources) == 0) return(list())

  out <- lapply(enabled_sources, function(src) {
    dat <- load_single_pqtl_source(src, target_genes = target_genes)
    list(
      source_name = src$source_name,
      source_id = src$source_id %||% src$source_name,
      sample_size = src$sample_size %||% NA_real_,
      ancestry = src$ancestry %||% NA_character_,
      source_role = src$source_role %||% NA_character_,
      protocol_status = src$protocol_status %||% NA_character_,
      outcome_scope = unlist(src$outcome_scope %||% character(0), use.names = FALSE),
      data = dat
    )
  })

  out
}

make_pqtl_source_summary <- function(pqtl_sources) {
  if (length(pqtl_sources) == 0) return(data.frame())

  bind_rows(lapply(pqtl_sources, function(src) {
    dat <- src$data
    data.frame(
      pqtl_source = src$source_name,
      pqtl_source_id = src$source_id,
      sample_size = src$sample_size,
      ancestry = src$ancestry,
      source_role = src$source_role %||% NA_character_,
      protocol_status = src$protocol_status %||% NA_character_,
      outcome_scope = if (length(src$outcome_scope %||% character(0)) == 0) {
        NA_character_
      } else {
        paste(src$outcome_scope, collapse = ";")
      },
      n_rows = if (is.null(dat)) 0L else nrow(dat),
      n_genes = if (is.null(dat)) 0L else dplyr::n_distinct(dat$gene),
      n_snps = if (is.null(dat)) 0L else dplyr::n_distinct(dat$SNP),
      stringsAsFactors = FALSE
    )
  }))
}
