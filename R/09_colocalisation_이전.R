# R/09_colocalisation.R
# Colocalisation analysis utilities
# Version: 1.0
# Author: Mihye Kwon

suppressPackageStartupMessages({
    library(dplyr)
    library(TwoSampleMR)
    library(coloc)
})

`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
}


normalise_localisation_config <- function(protocol_cfg, run_cfg = NULL) {
    if (!is.null(protocol_cfg)) {
        if (is.null(protocol_cfg$colocalisation) && !is.null(protocol_cfg$colocalization)) {
            protocol_cfg$colocalisation <- protocol_cfg$colocalization
        }
        if (is.null(protocol_cfg$colocalization) && !is.null(protocol_cfg$colocalisation)) {
            protocol_cfg$colocalization <- protocol_cfg$colocalisation
        }
    }
    if (!is.null(run_cfg)) {
        if (is.null(run_cfg$colocalisation_run) && !is.null(run_cfg$colocalization_run)) {
            run_cfg$colocalisation_run <- run_cfg$colocalization_run
        }
        if (is.null(run_cfg$colocalization_run) && !is.null(run_cfg$colocalisation_run)) {
            run_cfg$colocalization_run <- run_cfg$colocalisation_run
        }
    }
    list(protocol_cfg = protocol_cfg, run_cfg = run_cfg)
}


get_coloc_outcome_objects <- function(run_cfg, outcomes) {
    mode <- run_cfg$execution_scope$colocalisation_mode %||%
        "priority_targets_primary_outcome"
    
    if (identical(mode, "priority_targets_primary_outcome")) {
        return(Filter(function(x) identical(x$role, "primary"), outcomes))
    }
    
    if (identical(mode, "priority_targets_multi_outcome")) {
        core_names <- run_cfg$colocalisation_run$coloc_outcomes$core %||% character(0)
        supp_names <- run_cfg$colocalisation_run$coloc_outcomes$supplementary %||% character(0)
        target_names <- c(core_names, supp_names)
        matched <- Filter(function(x) x$name %in% target_names, outcomes)
        if (length(matched) == 0) {
            # fallback: primary만
            return(Filter(function(x) identical(x$role, "primary"), outcomes))
        }
        return(matched)
    }
    
    # 기타 알 수 없는 mode: fallback primary
    Filter(function(x) identical(x$role, "primary"), outcomes)
}

safe_coloc_log <- function(file = NULL, message) {
    if (is.null(file) || !nzchar(file)) return(invisible(NULL))
    if (exists("log_message", mode = "function", inherits = TRUE)) {
        log_message(file, message)
    }
    invisible(NULL)
}

classify_coloc_support <- function(pp_h3, pp_h4, h4_cutoff = 0.8, h3_cutoff = 0.5) {
    if (is.na(pp_h4)) return("Not_run")
    if (pp_h4 >= h4_cutoff) return("Strong_coloc")
    if (!is.na(pp_h3) && pp_h3 > pp_h4 && pp_h3 >= h3_cutoff) return("Distinct_signals_likely")
    "Weak_or_uncertain"
}

collapse_notes <- function(...) {
    x <- unlist(list(...), use.names = FALSE)
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0) return(NA_character_)
    paste(x, collapse = " | ")
}

empty_gwas_region_for_coloc <- function() {
    data.frame(
        snp = character(0),
        beta_gwas = numeric(0),
        se_gwas = numeric(0),
        varbeta_gwas = numeric(0),
        pval_gwas = numeric(0),
        eaf_gwas = numeric(0),
        effect_allele_outcome = character(0),
        other_allele_outcome = character(0),
        samplesize_gwas = numeric(0),
        ncase_gwas = numeric(0),
        ncontrol_gwas = numeric(0),
        stringsAsFactors = FALSE
    )
}

attach_fetch_status <- function(df, status) {
    attr(df, "coloc_fetch_status") <- status
    df
}

make_fetch_status <- function(requested_snps = 0L,
    chunk_size = NA_integer_,
    total_chunks = 0L,
    successful_chunks = 0L,
    failed_chunks = 0L,
    raw_rows = 0L,
    usable_rows = 0L,
    reason = NA_character_,
    note = NA_character_,
    failed_chunk_messages = character(0)) {
    list(
        requested_snps = as.integer(requested_snps),
        chunk_size = as.integer(chunk_size),
        total_chunks = as.integer(total_chunks),
        successful_chunks = as.integer(successful_chunks),
        failed_chunks = as.integer(failed_chunks),
        raw_rows = as.integer(raw_rows),
        usable_rows = as.integer(usable_rows),
        reason = reason,
        note = note,
        failed_chunk_messages = failed_chunk_messages
    )
}

make_coloc_result_row <- function(gene,
    outcome_obj,
    eqtl_region,
    nsnps_overlap = 0,
    n_gwas_representative = NA_real_,
    PP.H0 = NA_real_,
    PP.H1 = NA_real_,
    PP.H2 = NA_real_,
    PP.H3 = NA_real_,
    PP.H4 = NA_real_,
    coloc_support = "Not_run",
    note = NA_character_,
    first_pass_method = "coloc_abf",
    method = "coloc_abf",
    follow_up_method = NA_character_,
    abf_PP.H0 = NA_real_,
    abf_PP.H1 = NA_real_,
    abf_PP.H2 = NA_real_,
    abf_PP.H3 = NA_real_,
    abf_PP.H4 = NA_real_,
    abf_coloc_support = NA_character_,
    final_PP.H4 = PP.H4,
    final_coloc_support = coloc_support,
    susie_triggered = FALSE,
    susie_trigger_reason = NA_character_,
    susie_attempted = FALSE,
    susie_status = "not_attempted",
    susie_signal_pairs = NA_integer_,
    susie_PP.H3_max = NA_real_,
    susie_PP.H4_max = NA_real_,
    ld_population = NA_character_,
    ld_reference_panel = NA_character_,
    ld_source = NA_character_,
    ld_matrix_status = NA_character_,
    ld_matrix_note = NA_character_) {
    data.frame(
        gene = gene,
        outcome_name = outcome_obj$name,
        outcome_id = outcome_obj$id,
        region_chr = if (nrow(eqtl_region) > 0) unique(eqtl_region$chr)[1] else NA,
        region_start = if (nrow(eqtl_region) > 0) min(eqtl_region$pos, na.rm = TRUE) else NA,
        region_end = if (nrow(eqtl_region) > 0) max(eqtl_region$pos, na.rm = TRUE) else NA,
        first_pass_method = first_pass_method,
        follow_up_method = follow_up_method,
        method = method,
        nsnps_eqtl = nrow(eqtl_region),
        nsnps_overlap = nsnps_overlap,
        n_eqtl_representative = if (nrow(eqtl_region) > 0) unique(eqtl_region$n_eqtl_gene_median)[1] else NA_real_,
        n_gwas_representative = n_gwas_representative,
        PP.H0 = PP.H0,
        PP.H1 = PP.H1,
        PP.H2 = PP.H2,
        PP.H3 = PP.H3,
        PP.H4 = PP.H4,
        abf_PP.H0 = abf_PP.H0,
        abf_PP.H1 = abf_PP.H1,
        abf_PP.H2 = abf_PP.H2,
        abf_PP.H3 = abf_PP.H3,
        abf_PP.H4 = abf_PP.H4,
        abf_coloc_support = abf_coloc_support,
        final_PP.H4 = final_PP.H4,
        coloc_support = coloc_support,
        final_coloc_support = final_coloc_support,
        susie_triggered = susie_triggered,
        susie_trigger_reason = susie_trigger_reason,
        susie_attempted = susie_attempted,
        susie_status = susie_status,
        susie_signal_pairs = susie_signal_pairs,
        susie_PP.H3_max = susie_PP.H3_max,
        susie_PP.H4_max = susie_PP.H4_max,
        ld_population = ld_population,
        ld_reference_panel = ld_reference_panel,
        ld_source = ld_source,
        ld_matrix_status = ld_matrix_status,
        ld_matrix_note = ld_matrix_note,
        note = note,
        stringsAsFactors = FALSE
    )
}

prepare_eqtl_region_for_coloc <- function(exp_raw, target_gene) {
    dat <- exp_raw %>%
        dplyr::filter(.data$gene == tolower(.env$target_gene)) %>%
        dplyr::transmute(
            gene = .data$gene,
            snp = .data$SNP,
            chr = .data$chr,
            pos = .data$pos,
            beta_eqtl = .data$beta,
            se_eqtl = .data$se,
            varbeta_eqtl = .data$se^2,
            pval_eqtl = .data$pval,
            maf = .data$maf,
            n_eqtl_variant = .data$n_eqtl_variant,
            n_eqtl_gene_median = .data$n_eqtl_gene_median,
            n_eqtl_gene_max = .data$n_eqtl_gene_max,
            type_eqtl = "quant"
        ) %>%
        dplyr::filter(
            !is.na(.data$snp),
            nzchar(.data$snp),
            !is.na(.data$chr),
            !is.na(.data$pos),
            !is.na(.data$beta_eqtl),
            !is.na(.data$varbeta_eqtl),
            !is.na(.data$pval_eqtl),
            !is.na(.data$maf)
        ) %>%
        dplyr::arrange(.data$pval_eqtl, dplyr::desc(.data$n_eqtl_variant)) %>%
        dplyr::distinct(.data$snp, .keep_all = TRUE)
    
    dat
}

clean_outcome_region_for_coloc <- function(out_dat) {
    if (is.null(out_dat) || nrow(out_dat) == 0) {
        return(empty_gwas_region_for_coloc())
    }
    
    samplesize_vec <- if ("samplesize.outcome" %in% names(out_dat)) {
        suppressWarnings(as.numeric(out_dat$samplesize.outcome))
    } else {
        rep(NA_real_, nrow(out_dat))
    }
    
    ncase_vec <- if ("ncase.outcome" %in% names(out_dat)) {
        suppressWarnings(as.numeric(out_dat$ncase.outcome))
    } else {
        rep(NA_real_, nrow(out_dat))
    }
    
    ncontrol_vec <- if ("ncontrol.outcome" %in% names(out_dat)) {
        suppressWarnings(as.numeric(out_dat$ncontrol.outcome))
    } else {
        rep(NA_real_, nrow(out_dat))
    }
    
    out_dat %>%
        dplyr::mutate(
            `..samplesize_gwas` = samplesize_vec,
            `..ncase_gwas` = ncase_vec,
            `..ncontrol_gwas` = ncontrol_vec
        ) %>%
        dplyr::transmute(
            snp = .data$SNP,
            beta_gwas = .data$beta.outcome,
            se_gwas = .data$se.outcome,
            varbeta_gwas = .data$se.outcome^2,
            pval_gwas = .data$pval.outcome,
            eaf_gwas = .data$eaf.outcome,
            effect_allele_outcome = .data$effect_allele.outcome,
            other_allele_outcome = .data$other_allele.outcome,
            samplesize_gwas = .data$`..samplesize_gwas`,
            ncase_gwas = .data$`..ncase_gwas`,
            ncontrol_gwas = .data$`..ncontrol_gwas`
        ) %>%
        dplyr::filter(
            !is.na(.data$snp),
            nzchar(.data$snp),
            !is.na(.data$beta_gwas),
            !is.na(.data$varbeta_gwas),
            !is.na(.data$pval_gwas)
        ) %>%
        dplyr::arrange(.data$pval_gwas) %>%
        dplyr::distinct(.data$snp, .keep_all = TRUE)
}

get_outcome_cc_meta <- function(gwas_region, outcome_obj, protocol_cfg = NULL) {
    pick_first_numeric <- function(x) {
        x <- suppressWarnings(as.numeric(x))
        x <- x[is.finite(x) & x > 0]
        if (length(x) == 0) return(NA_real_)
        x[1]
    }
    
    as_outcome_list <- function(x) {
        if (is.null(x)) return(list())
        if (is.list(x) && !is.null(x$id)) return(list(x))
        if (is.list(x)) return(unname(x))
        list()
    }
    
    n_cases <- pick_first_numeric(outcome_obj$n_cases)
    n_controls <- pick_first_numeric(outcome_obj$n_controls)
    
    if ((!is.finite(n_cases) || !is.finite(n_controls)) && !is.null(gwas_region) && nrow(gwas_region) > 0) {
        if ("ncase_gwas" %in% names(gwas_region)) {
            n_cases <- pick_first_numeric(gwas_region$ncase_gwas)
        }
        if ("ncontrol_gwas" %in% names(gwas_region)) {
            n_controls <- pick_first_numeric(gwas_region$ncontrol_gwas)
        }
    }
    
    if ((!is.finite(n_cases) || !is.finite(n_controls)) && !is.null(protocol_cfg)) {
        all_outcomes <- c(
            as_outcome_list(protocol_cfg$outcomes$primary %||% NULL),
            as_outcome_list(protocol_cfg$outcomes$replication %||% NULL),
            as_outcome_list(protocol_cfg$outcomes$secondary %||% NULL)
        )
        
        hit <- Filter(
            function(x) identical(x$id %||% NA_character_, outcome_obj$id %||% NA_character_),
            all_outcomes
        )
        
        if (length(hit) >= 1) {
            n_cases <- pick_first_numeric(hit[[1]]$n_cases)
            n_controls <- pick_first_numeric(hit[[1]]$n_controls)
        }
    }
    
    if (!is.finite(n_cases) || !is.finite(n_controls) || n_cases <= 0 || n_controls <= 0) {
        stop(paste0("Missing numeric n_cases/n_controls for outcome ", outcome_obj$id %||% NA_character_))
    }
    
    list(
        n_cases = n_cases,
        n_controls = n_controls,
        N = n_cases + n_controls,
        s = n_cases / (n_cases + n_controls)
    )
}

extract_outcome_chunk_exact <- function(snps, outcome_id, err_file = NULL) {
    tryCatch(
        TwoSampleMR::extract_outcome_data(
            snps = snps,
            outcomes = outcome_id,
            proxies = FALSE
        ),
        error = function(e) {
            safe_coloc_log(
                err_file,
                paste0(
                    "COLOC OUTCOME EXTRACT CHUNK FAIL | outcome_id=", outcome_id,
                    " requested_snps=", length(snps),
                    " msg=", e$message
                )
            )
            structure(list(message = e$message), class = "coloc_extract_error")
        }
    )
}

extract_outcome_region_for_coloc <- function(eqtl_region,
    outcome_id,
    err_file = NULL,
    chunk_size = 100L,
    max_retries = 3L,
    retry_wait_sec = 1) {
    snps <- unique(eqtl_region$snp)
    snps <- snps[!is.na(snps) & nzchar(snps)]
    
    if (length(snps) == 0) {
        status <- make_fetch_status(
            requested_snps = 0L,
            chunk_size = chunk_size,
            total_chunks = 0L,
            reason = "no_requested_snps",
            note = "No valid requested SNPs for exact-only outcome extraction"
        )
        return(attach_fetch_status(empty_gwas_region_for_coloc(), status))
    }
    
    snp_chunks <- split(snps, ceiling(seq_along(snps) / chunk_size))
    total_chunks <- length(snp_chunks)
    raw_chunks <- vector("list", total_chunks)
    failed_messages <- character(0)
    successful_chunks <- 0L
    failed_chunks <- 0L
    
    safe_coloc_log(
        err_file,
        paste0(
            "COLOC OUTCOME EXTRACT START | outcome_id=", outcome_id,
            " requested_snps=", length(snps),
            " chunk_size=", chunk_size,
            " total_chunks=", total_chunks,
            " mode=exact_only"
        )
    )
    
    for (i in seq_along(snp_chunks)) {
        chunk_snps <- snp_chunks[[i]]
        chunk_result <- NULL
        chunk_error_message <- NULL
        
        for (attempt in seq_len(max_retries)) {
            safe_coloc_log(
                err_file,
                paste0(
                    "COLOC OUTCOME EXTRACT ATTEMPT | outcome_id=", outcome_id,
                    " chunk=", i, "/", total_chunks,
                    " attempt=", attempt,
                    " requested_snps=", length(chunk_snps)
                )
            )
            
            chunk_result <- extract_outcome_chunk_exact(
                snps = chunk_snps,
                outcome_id = outcome_id,
                err_file = err_file
            )
            
            if (!inherits(chunk_result, "coloc_extract_error")) {
                break
            }
            
            chunk_error_message <- chunk_result$message
            if (attempt < max_retries) {
                Sys.sleep(retry_wait_sec)
            }
        }
        
        if (inherits(chunk_result, "coloc_extract_error")) {
            failed_chunks <- failed_chunks + 1L
            failed_messages <- c(
                failed_messages,
                paste0("chunk ", i, "/", total_chunks, ": ", chunk_error_message)
            )
            next
        }
        
        successful_chunks <- successful_chunks + 1L
        if (!is.null(chunk_result) && nrow(chunk_result) > 0) {
            raw_chunks[[i]] <- chunk_result
        }
    }
    
    raw_bound <- dplyr::bind_rows(raw_chunks)
    cleaned <- clean_outcome_region_for_coloc(raw_bound)
    
    reason <- NA_character_
    note <- NA_character_
    
    if (successful_chunks == 0L && failed_chunks > 0L) {
        reason <- "outcome_query_failed_all_chunks"
        note <- paste0(
            "Exact-only outcome extraction failed in all ", total_chunks,
            " chunks"
        )
    } else if (nrow(raw_bound) == 0L && failed_chunks == 0L) {
        reason <- "outcome_query_zero_rows_no_error"
        note <- "Exact-only outcome extraction returned zero rows"
    } else if (nrow(cleaned) == 0L && nrow(raw_bound) > 0L) {
        reason <- "outcome_rows_filtered_to_zero"
        note <- paste0(
            "Outcome extraction returned rows but zero usable rows remained after cleaning",
            if (failed_chunks > 0L) paste0("; ", failed_chunks, "/", total_chunks, " chunks failed") else ""
        )
    } else if (failed_chunks > 0L) {
        reason <- "outcome_query_partial_failure_but_rows_recovered"
        note <- paste0(
            "Exact-only retrieval with partial recovery; ",
            failed_chunks, "/", total_chunks, " chunks failed"
        )
    } else {
        reason <- "exact_only_success"
        note <- "Exact-only outcome extraction succeeded"
    }
    
    if (length(failed_messages) > 0L) {
        safe_coloc_log(
            err_file,
            paste0(
                "COLOC OUTCOME EXTRACT FAILURES | outcome_id=", outcome_id,
                " details=", paste(failed_messages, collapse = " || ")
            )
        )
    }
    
    safe_coloc_log(
        err_file,
        paste0(
            "COLOC OUTCOME EXTRACT END | outcome_id=", outcome_id,
            " requested_snps=", length(snps),
            " raw_rows=", nrow(raw_bound),
            " usable_rows=", nrow(cleaned),
            " successful_chunks=", successful_chunks,
            " failed_chunks=", failed_chunks,
            " reason=", reason
        )
    )
    
    status <- make_fetch_status(
        requested_snps = length(snps),
        chunk_size = chunk_size,
        total_chunks = total_chunks,
        successful_chunks = successful_chunks,
        failed_chunks = failed_chunks,
        raw_rows = nrow(raw_bound),
        usable_rows = nrow(cleaned),
        reason = reason,
        note = note,
        failed_chunk_messages = failed_messages
    )
    
    attach_fetch_status(cleaned, status)
}


get_ld_population_for_outcome <- function(outcome_obj, protocol_cfg) {
    anc <- toupper(outcome_obj$ancestry %||% outcome_obj$population %||% "EUR")
    if (identical(anc, "EAS")) return("EAS")
    "EUR"
}

get_ld_reference_panel_name <- function(protocol_cfg, ld_population) {
    protocol_cfg$instrument_selection$ld_reference[[ld_population]]$panel %||%
        protocol_cfg$instrument_selection$ld_reference[[ld_population]] %||%
        paste0("1000G_phase3_", ld_population)
}

get_ld_local_bfile <- function(protocol_cfg, ld_population) {
    protocol_cfg$colocalisation$ld_matrix$local_plink$reference_bfiles[[ld_population]] %||% ""
}

get_ld_local_plink_bin <- function(protocol_cfg) {
    protocol_cfg$colocalisation$ld_matrix$local_plink$plink_bin %||% ""
}


get_case_prop_for_outcome <- function(outcome_obj) {
    n_cases <- suppressWarnings(as.numeric(outcome_obj$n_cases %||% NA_real_))
    n_controls <- suppressWarnings(as.numeric(outcome_obj$n_controls %||% NA_real_))
    if (!is.finite(n_cases) || !is.finite(n_controls) || (n_cases + n_controls) <= 0) {
        return(NA_real_)
    }
    n_cases / (n_cases + n_controls)
}

normalise_ld_matrix <- function(ld) {
    if (is.null(ld)) return(NULL)
    ld <- as.matrix(ld)
    storage.mode(ld) <- "numeric"
    rn <- rownames(ld)
    cn <- colnames(ld)
    if (is.null(rn) || is.null(cn)) return(NULL)
    keep <- intersect(rn, cn)
    if (length(keep) < 2) return(NULL)
    ld[keep, keep, drop = FALSE]
}

compute_ld_matrix_for_coloc <- function(snps, ld_population, protocol_cfg, err_file = NULL) {
    snps <- unique(snps)
    snps <- snps[!is.na(snps) & nzchar(snps)]
    
    if (length(snps) < 2) {
        return(list(
            ok = FALSE,
            status = "insufficient_snps",
            note = "Need at least 2 SNPs to compute local LD matrix",
            source = "local_plink",
            ld_population = ld_population,
            reference_panel = get_ld_reference_panel_name(protocol_cfg, ld_population),
            ld = NULL
        ))
    }
    
    local_enabled <- isTRUE(protocol_cfg$colocalisation$ld_matrix$local_plink$enabled)
    bfile <- get_ld_local_bfile(protocol_cfg, ld_population)
    plink_bin <- get_ld_local_plink_bin(protocol_cfg)
    
    if (!local_enabled) {
        return(list(
            ok = FALSE,
            status = "local_ld_disabled",
            note = "Local LD is required for coloc.susie but local_plink is disabled",
            source = "local_plink",
            ld_population = ld_population,
            reference_panel = get_ld_reference_panel_name(protocol_cfg, ld_population),
            ld = NULL
        ))
    }
    
    if (!nzchar(bfile) || !file.exists(paste0(bfile, ".bed")) || !file.exists(paste0(bfile, ".bim")) || !file.exists(paste0(bfile, ".fam"))) {
        return(list(
            ok = FALSE,
            status = "missing_local_reference",
            note = "Local LD reference bfile is missing (.bed/.bim/.fam)",
            source = "local_plink",
            ld_population = ld_population,
            reference_panel = get_ld_reference_panel_name(protocol_cfg, ld_population),
            ld = NULL
        ))
    }
    
    if (!nzchar(plink_bin) || !file.exists(plink_bin)) {
        return(list(
            ok = FALSE,
            status = "missing_plink_binary",
            note = "PLINK binary not found for local LD computation",
            source = "local_plink",
            ld_population = ld_population,
            reference_panel = get_ld_reference_panel_name(protocol_cfg, ld_population),
            ld = NULL
        ))
    }
    
    if (!requireNamespace("ieugwasr", quietly = TRUE)) {
        return(list(
            ok = FALSE,
            status = "missing_ieugwasr_package",
            note = "ieugwasr package is required for ld_matrix_local",
            source = "local_plink",
            ld_population = ld_population,
            reference_panel = get_ld_reference_panel_name(protocol_cfg, ld_population),
            ld = NULL
        ))
    }
    
    ld <- tryCatch(
        ieugwasr::ld_matrix_local(variants = snps, bfile = bfile, plink_bin = plink_bin, with_alleles = FALSE),
        error = function(e) e
    )
    
    if (inherits(ld, "error")) {
        return(list(
            ok = FALSE,
            status = "local_ld_failed",
            note = paste("Local LD computation failed:", ld$message),
            source = "local_plink",
            ld_population = ld_population,
            reference_panel = get_ld_reference_panel_name(protocol_cfg, ld_population),
            ld = NULL
        ))
    }
    
    ld <- normalise_ld_matrix(ld)
    if (is.null(ld) || nrow(ld) < 2) {
        return(list(
            ok = FALSE,
            status = "local_ld_empty",
            note = "Local LD returned fewer than 2 SNPs after normalization",
            source = "local_plink",
            ld_population = ld_population,
            reference_panel = get_ld_reference_panel_name(protocol_cfg, ld_population),
            ld = NULL
        ))
    }
    
    list(
        ok = TRUE,
        status = "ok",
        note = "LD matrix obtained via local PLINK reference",
        source = "local_plink",
        ld_population = ld_population,
        reference_panel = get_ld_reference_panel_name(protocol_cfg, ld_population),
        ld = ld
    )
}


gene_in_cluster_or_shared_locus <- function(gene, protocol_cfg) {
    gene_lower <- tolower(gene %||% "")
    
    cluster_groups <- protocol_cfg$cluster_groups %||% list()
    cluster_hit_names <- vapply(
        Filter(function(x) gene_lower %in% tolower(x$genes %||% character(0)), cluster_groups),
        function(x) x$name %||% NA_character_,
        character(1)
    )
    
    locus_rules <- protocol_cfg$locus_specific_rules %||% list()
    shared_hit_names <- vapply(
        Filter(function(x) gene_lower %in% tolower(x$genes %||% character(0)), locus_rules),
        function(x) x$name %||% NA_character_,
        character(1)
    )
    
    list(
        cluster_hit = length(cluster_hit_names) > 0,
        shared_hit = length(shared_hit_names) > 0,
        cluster_names = unname(cluster_hit_names),
        shared_rule_names = unname(shared_hit_names)
    )
}

is_primary_mr_supportive <- function(gene, preferred_df, run_cfg = NULL) {
    if (is.null(preferred_df) || !is.data.frame(preferred_df) || nrow(preferred_df) == 0) {
        return(FALSE)
    }
    
    nm <- names(preferred_df)
    gene_col <- intersect(c("gene", "target", "target_gene", "exposure", "exposure_gene"), nm)[1]
    if (is.na(gene_col)) return(FALSE)
    
    d <- preferred_df[tolower(as.character(preferred_df[[gene_col]])) == tolower(gene), , drop = FALSE]
    if (nrow(d) == 0) return(FALSE)
    
    role_col <- intersect(c("outcome_role", "role", "analysis_role", "analysis_set", "outcome_set"), nm)[1]
    if (!is.na(role_col)) {
        primary_patterns <- c("primary", "eur_primary", "ra_eur_primary")
        keep <- grepl(paste(primary_patterns, collapse = "|"),
            tolower(as.character(d[[role_col]])))
        if (any(keep, na.rm = TRUE)) d <- d[keep, , drop = FALSE]
    }
    
    class_col <- intersect(c("benchmark_class", "target_class", "classification",
        "support_class", "mr_support", "support"), nm)[1]
    if (!is.na(class_col)) {
        supportive_class <- c("strong", "weak_or_partial", "supported", "suggestive")
        vals <- tolower(as.character(d[[class_col]]))
        if (any(vals %in% supportive_class, na.rm = TRUE)) return(TRUE)
    }
    
    q_col <- intersect(c("qval", "q_value", "fdr", "fdr_q"), nm)[1]
    if (!is.na(q_col)) {
        qv <- suppressWarnings(as.numeric(d[[q_col]]))
        if (any(is.finite(qv) & qv <= 0.10, na.rm = TRUE)) return(TRUE)
    }
    
    p_col <- intersect(c("pval", "p_value", "p", "mr_pval"), nm)[1]
    if (!is.na(p_col)) {
        pv <- suppressWarnings(as.numeric(d[[p_col]]))
        if (any(is.finite(pv) & pv <= 0.05, na.rm = TRUE)) return(TRUE)
    }
    
    FALSE
}

should_run_coloc_susie <- function(gene, abf_support, preferred_df, run_cfg, protocol_cfg) {
    if (!isTRUE(run_cfg$colocalisation_run$allow_extended_method_today %||% FALSE)) {
        return(list(run = FALSE, reason = "run_profile_disables_extended_coloc"))
    }
    
    if (!identical(protocol_cfg$colocalisation$follow_up_method %||% "", "coloc_susie")) {
        return(list(run = FALSE, reason = "protocol_follow_up_method_not_susie"))
    }
    
    flags <- gene_in_cluster_or_shared_locus(gene, protocol_cfg)
    complex_ld_flag <- isTRUE(protocol_cfg$colocalisation$region_definition$complex_ld_region_caution) &&
        (flags$cluster_hit || flags$shared_hit)
    mr_supportive <- is_primary_mr_supportive(gene, preferred_df, run_cfg)
    
    reasons <- character(0)
    if (flags$cluster_hit || flags$shared_hit) reasons <- c(reasons, "shared_locus_or_cluster_region")
    if (complex_ld_flag) reasons <- c(reasons, "complex_ld_region")
    if (identical(abf_support, "Weak_or_uncertain") && mr_supportive) reasons <- c(reasons, "abf_inconclusive_but_mr_supportive")
    if (identical(abf_support, "Distinct_signals_likely")) reasons <- c(reasons, "abf_distinct_signals_likely")
    
    if (length(reasons) == 0) {
        return(list(run = FALSE, reason = "no_susie_trigger_met"))
    }
    
    list(run = TRUE, reason = paste(unique(reasons), collapse = "; "))
}

summarise_susie_result <- function(susie_res, protocol_cfg) {
    if (is.null(susie_res) || is.null(susie_res$summary) || nrow(as.data.frame(susie_res$summary)) == 0) {
        return(list(
            signal_pairs = 0L,
            PP.H0 = NA_real_,
            PP.H1 = NA_real_,
            PP.H2 = NA_real_,
            PP.H3 = NA_real_,
            PP.H4 = NA_real_,
            coloc_support = "Not_run"
        ))
    }
    
    summ <- as.data.frame(susie_res$summary)
    pp_h3 <- if ("PP.H3.abf" %in% names(summ)) max(summ$PP.H3.abf, na.rm = TRUE) else NA_real_
    pp_h4 <- if ("PP.H4.abf" %in% names(summ)) max(summ$PP.H4.abf, na.rm = TRUE) else NA_real_
    
    list(
        signal_pairs = nrow(summ),
        PP.H0 = if ("PP.H0.abf" %in% names(summ)) max(summ$PP.H0.abf, na.rm = TRUE) else NA_real_,
        PP.H1 = if ("PP.H1.abf" %in% names(summ)) max(summ$PP.H1.abf, na.rm = TRUE) else NA_real_,
        PP.H2 = if ("PP.H2.abf" %in% names(summ)) max(summ$PP.H2.abf, na.rm = TRUE) else NA_real_,
        PP.H3 = pp_h3,
        PP.H4 = pp_h4,
        coloc_support = classify_coloc_support(
            pp_h3 = pp_h3,
            pp_h4 = pp_h4,
            h4_cutoff = protocol_cfg$colocalisation$strong_support_threshold$pp_h4_ge
        )
    )
}


run_single_target_coloc_abf <- function(eqtl_region,
    gwas_region,
    gene,
    outcome_obj,
    protocol_cfg,
    run_cfg,
    preferred_df = NULL,
    fetch_status = NULL,
    err_file = NULL) {
    cfg_norm <- normalise_localisation_config(protocol_cfg, run_cfg)
    protocol_cfg <- cfg_norm$protocol_cfg
    run_cfg <- cfg_norm$run_cfg
    merged <- eqtl_region %>%
        dplyr::inner_join(gwas_region, by = "snp") %>%
        dplyr::filter(
            !is.na(.data$snp),
            nzchar(.data$snp),
            !is.na(.data$beta_eqtl),
            !is.na(.data$varbeta_eqtl),
            !is.na(.data$beta_gwas),
            !is.na(.data$varbeta_gwas),
            !is.na(.data$maf)
        ) %>%
        dplyr::arrange(.data$pval_eqtl, .data$pval_gwas) %>%
        dplyr::distinct(.data$snp, .keep_all = TRUE)
    
    n_gwas_rep <- suppressWarnings(max(merged$samplesize_gwas, na.rm = TRUE))
    if (!is.finite(n_gwas_rep)) n_gwas_rep <- NA_real_
    
    empty_row <- function(note_text, nsnps_overlap = nrow(merged), susie_triggered = FALSE,
        susie_trigger_reason = NA_character_, ld_info = NULL,
        susie_attempted = FALSE, susie_status = "not_attempted") {
        make_coloc_result_row(
            gene = gene,
            outcome_obj = outcome_obj,
            eqtl_region = eqtl_region,
            nsnps_overlap = nsnps_overlap,
            n_gwas_representative = n_gwas_rep,
            coloc_support = "Not_run",
            final_coloc_support = "Not_run",
            note = collapse_notes(
                note_text,
                if (!is.null(fetch_status)) fetch_status$note else NA_character_
            ),
            susie_triggered = susie_triggered,
            susie_trigger_reason = susie_trigger_reason,
            susie_attempted = susie_attempted,
            susie_status = susie_status,
            ld_population = ld_info$ld_population %||% NA_character_,
            ld_reference_panel = ld_info$reference_panel %||% NA_character_,
            ld_source = ld_info$source %||% NA_character_,
            ld_matrix_status = ld_info$status %||% NA_character_,
            ld_matrix_note = ld_info$note %||% NA_character_
        )
    }
    
    if (anyDuplicated(merged$snp)) {
        return(empty_row("Duplicated SNPs remain after deduplication"))
    }
    
    if (nrow(merged) == 0) {
        return(empty_row("No exact overlapping SNPs after merge/cleaning", nsnps_overlap = 0))
    }
    
    n_eqtl_rep <- unique(merged$n_eqtl_gene_median)[1]
    if (!is.finite(n_eqtl_rep)) n_eqtl_rep <- NA_real_
    
    cc_meta <- get_outcome_cc_meta(gwas_region, outcome_obj, protocol_cfg)
    n_gwas_rep <- cc_meta$N
    
    d1 <- list(
        beta = merged$beta_eqtl,
        varbeta = merged$varbeta_eqtl,
        snp = merged$snp,
        position = merged$pos,
        type = "quant",
        N = n_eqtl_rep,
        MAF = merged$maf
    )
    
    d2 <- list(
        beta = merged$beta_gwas,
        varbeta = merged$varbeta_gwas,
        snp = merged$snp,
        position = merged$pos,
        type = "cc",
        N = cc_meta$N,
        s = cc_meta$s,
        MAF = merged$maf
    )
    
    abf_res <- tryCatch(
        coloc::coloc.abf(
            dataset1 = d1,
            dataset2 = d2,
            p1 = protocol_cfg$colocalisation$priors$p1,
            p2 = protocol_cfg$colocalisation$priors$p2,
            p12 = protocol_cfg$colocalisation$priors$p12
        ),
        error = function(e) e
    )
    
    if (inherits(abf_res, "error")) {
        return(empty_row(paste("coloc.abf failed:", abf_res$message)))
    }
    
    abf_summ <- abf_res$summary
    abf_pp_h3 <- unname(abf_summ["PP.H3.abf"])
    abf_pp_h4 <- unname(abf_summ["PP.H4.abf"])
    abf_support <- classify_coloc_support(
        pp_h3 = abf_pp_h3,
        pp_h4 = abf_pp_h4,
        h4_cutoff = protocol_cfg$colocalisation$strong_support_threshold$pp_h4_ge
    )
    
    caution_note <- NA_character_
    if (tolower(gene) %in% c("tnf", "lta")) {
        caution_note <- "Complex shared locus / MHC-proximal caution"
    }
    
    final_vals <- list(
        method = "coloc_abf",
        PP.H0 = unname(abf_summ["PP.H0.abf"]),
        PP.H1 = unname(abf_summ["PP.H1.abf"]),
        PP.H2 = unname(abf_summ["PP.H2.abf"]),
        PP.H3 = abf_pp_h3,
        PP.H4 = abf_pp_h4,
        coloc_support = abf_support,
        final_PP.H4 = abf_pp_h4,
        final_coloc_support = abf_support,
        susie_triggered = FALSE,
        susie_trigger_reason = NA_character_,
        susie_attempted = FALSE,
        susie_status = "not_attempted",
        susie_signal_pairs = NA_integer_,
        susie_PP.H3_max = NA_real_,
        susie_PP.H4_max = NA_real_,
        ld_population = NA_character_,
        ld_reference_panel = NA_character_,
        ld_source = NA_character_,
        ld_matrix_status = NA_character_,
        ld_matrix_note = NA_character_,
        follow_up_method = NA_character_
    )
    
    susie_decision <- should_run_coloc_susie(
        gene = gene,
        abf_support = abf_support,
        preferred_df = preferred_df,
        run_cfg = run_cfg,
        protocol_cfg = protocol_cfg
    )
    final_vals$susie_triggered <- isTRUE(susie_decision$run)
    final_vals$susie_trigger_reason <- susie_decision$reason
    
    if (isTRUE(susie_decision$run)) {
        ld_population <- get_ld_population_for_outcome(outcome_obj, protocol_cfg)
        ld_info <- compute_ld_matrix_for_coloc(merged$snp, ld_population, protocol_cfg, err_file = err_file)
        
        final_vals$ld_population <- ld_info$ld_population %||% ld_population
        final_vals$ld_reference_panel <- ld_info$reference_panel %||% NA_character_
        final_vals$ld_source <- ld_info$source %||% NA_character_
        final_vals$ld_matrix_status <- ld_info$status %||% NA_character_
        final_vals$ld_matrix_note <- ld_info$note %||% NA_character_
        final_vals$susie_attempted <- TRUE
        final_vals$follow_up_method <- "coloc_susie"
        
        minimum_overlap <- as.integer(protocol_cfg$colocalisation$susie_follow_up$minimum_overlap_snps %||% 50L)
        
        if (nrow(merged) < minimum_overlap) {
            final_vals$susie_status <- "skipped_too_few_exact_overlap_snps"
            final_vals$ld_matrix_note <- collapse_notes(
                final_vals$ld_matrix_note,
                paste0("Only ", nrow(merged), " exact overlapping SNPs; minimum required for SuSiE is ", minimum_overlap)
            )
        } else if (is.null(ld_info$ld)) {
            final_vals$susie_status <- paste0("skipped_", ld_info$status %||% "ld_unavailable")
        } else {
            keep_snps <- intersect(merged$snp, rownames(ld_info$ld))
            keep_snps <- keep_snps[!is.na(keep_snps)]
            min_ld_snps <- as.integer(protocol_cfg$colocalisation$susie_follow_up$minimum_ld_snps %||% 20L)
            
            if (length(keep_snps) < min_ld_snps) {
                final_vals$susie_status <- "skipped_too_few_snps_after_ld_intersection"
                final_vals$ld_matrix_note <- collapse_notes(final_vals$ld_matrix_note,
                    paste0("Only ", length(keep_snps), " SNPs remained after LD intersection"))
            } else {
                merged_s <- merged %>%
                    dplyr::filter(.data$snp %in% keep_snps) %>%
                    dplyr::arrange(match(.data$snp, keep_snps))
                
                ld <- ld_info$ld[keep_snps, keep_snps, drop = FALSE]
                
                d1s <- list(
                    beta = merged_s$beta_eqtl,
                    varbeta = merged_s$varbeta_eqtl,
                    snp = merged_s$snp,
                    position = merged_s$pos,
                    type = "quant",
                    N = n_eqtl_rep,
                    MAF = merged_s$maf,
                    LD = ld
                )
                
                d2s <- list(
                    beta = merged_s$beta_gwas,
                    varbeta = merged_s$varbeta_gwas,
                    snp = merged_s$snp,
                    position = merged_s$pos,
                    type = "cc",
                    N = cc_meta$N,
                    s = cc_meta$s,
                    MAF = merged_s$maf,
                    LD = ld
                )
                
                susie_res <- tryCatch(
                    {
                        s1 <- coloc::runsusie(d1s)
                        s2 <- coloc::runsusie(d2s)
                        coloc::coloc.susie(
                            dataset1 = s1,
                            dataset2 = s2,
                            p1 = protocol_cfg$colocalisation$priors$p1,
                            p2 = protocol_cfg$colocalisation$priors$p2,
                            p12 = protocol_cfg$colocalisation$priors$p12
                        )
                    },
                    error = function(e) e
                )
                
                if (inherits(susie_res, "error")) {
                    final_vals$susie_status <- "failed"
                    final_vals$ld_matrix_note <- collapse_notes(final_vals$ld_matrix_note,
                        paste("coloc.susie failed:", susie_res$message))
                } else {
                    susie_sum <- summarise_susie_result(susie_res, protocol_cfg)
                    final_vals$method <- "coloc_susie"
                    final_vals$PP.H0 <- susie_sum$PP.H0
                    final_vals$PP.H1 <- susie_sum$PP.H1
                    final_vals$PP.H2 <- susie_sum$PP.H2
                    final_vals$PP.H3 <- susie_sum$PP.H3
                    final_vals$PP.H4 <- susie_sum$PP.H4
                    final_vals$coloc_support <- susie_sum$coloc_support
                    final_vals$final_PP.H4 <- susie_sum$PP.H4
                    final_vals$final_coloc_support <- susie_sum$coloc_support
                    final_vals$susie_signal_pairs <- susie_sum$signal_pairs
                    final_vals$susie_PP.H3_max <- susie_sum$PP.H3
                    final_vals$susie_PP.H4_max <- susie_sum$PP.H4
                    final_vals$susie_status <- "ok"
                }
            }
        }
    }
    
    make_coloc_result_row(
        gene = gene,
        outcome_obj = outcome_obj,
        eqtl_region = eqtl_region,
        nsnps_overlap = nrow(merged),
        n_gwas_representative = n_gwas_rep,
        PP.H0 = final_vals$PP.H0,
        PP.H1 = final_vals$PP.H1,
        PP.H2 = final_vals$PP.H2,
        PP.H3 = final_vals$PP.H3,
        PP.H4 = final_vals$PP.H4,
        coloc_support = final_vals$coloc_support,
        note = collapse_notes(
            caution_note,
            if (!is.null(fetch_status) && identical(fetch_status$reason, "outcome_query_partial_failure_but_rows_recovered")) {
                fetch_status$note
            } else {
                NA_character_
            },
            if (identical(final_vals$susie_status, "failed")) final_vals$ld_matrix_note else NA_character_
        ),
        method = final_vals$method,
        follow_up_method = final_vals$follow_up_method,
        abf_PP.H0 = unname(abf_summ["PP.H0.abf"]),
        abf_PP.H1 = unname(abf_summ["PP.H1.abf"]),
        abf_PP.H2 = unname(abf_summ["PP.H2.abf"]),
        abf_PP.H3 = abf_pp_h3,
        abf_PP.H4 = abf_pp_h4,
        abf_coloc_support = abf_support,
        final_PP.H4 = final_vals$final_PP.H4,
        final_coloc_support = final_vals$final_coloc_support,
        susie_triggered = final_vals$susie_triggered,
        susie_trigger_reason = final_vals$susie_trigger_reason,
        susie_attempted = final_vals$susie_attempted,
        susie_status = final_vals$susie_status,
        susie_signal_pairs = final_vals$susie_signal_pairs,
        susie_PP.H3_max = final_vals$susie_PP.H3_max,
        susie_PP.H4_max = final_vals$susie_PP.H4_max,
        ld_population = final_vals$ld_population,
        ld_reference_panel = final_vals$ld_reference_panel,
        ld_source = final_vals$ld_source,
        ld_matrix_status = final_vals$ld_matrix_status,
        ld_matrix_note = final_vals$ld_matrix_note
    )
}


run_priority_target_coloc_summary <- function(exp_raw,
    protocol_cfg,
    run_cfg,
    outcomes,
    preferred_df = NULL,
    err_file = NULL,
    log_file = NULL) {
    cfg_norm <- normalise_localisation_config(protocol_cfg, run_cfg)
    protocol_cfg <- cfg_norm$protocol_cfg
    run_cfg <- cfg_norm$run_cfg
    if (!isTRUE(run_cfg$execution_scope$run_colocalisation)) {
        return(data.frame())
    }
    
    priority_targets <- run_cfg$colocalisation_run$targets
    coloc_outcomes <- get_coloc_outcome_objects(run_cfg, outcomes)
    
    if (length(coloc_outcomes) == 0) {
        return(data.frame())
    }
    
    res_list <- list()
    
    for (outcome_obj in coloc_outcomes) {           # 외부 루프 열림
        safe_coloc_log(log_file, paste0("COLOC outcome loop start: ", outcome_obj$name))
        
        for (g in priority_targets) {                 # 내부 루프 열림
            safe_coloc_log(log_file, paste0("COLOC start gene=", g, " outcome=", outcome_obj$name))
            
            eqtl_region <- prepare_eqtl_region_for_coloc(exp_raw, g)
            
            if (nrow(eqtl_region) == 0) {
                res_list[[length(res_list) + 1L]] <- make_coloc_result_row(
                    gene = g,
                    outcome_obj = outcome_obj,
                    eqtl_region = eqtl_region,
                    coloc_support = "Not_run",
                    note = "No eQTL region data"
                )
                next
            }
            
            gwas_region <- extract_outcome_region_for_coloc(
                eqtl_region = eqtl_region,
                outcome_id = outcome_obj$id,
                err_file = err_file
            )
            fetch_status <- attr(gwas_region, "coloc_fetch_status")
            
            if (is.null(gwas_region) || nrow(gwas_region) == 0) {
                res_list[[length(res_list) + 1L]] <- make_coloc_result_row(
                    gene = g,
                    outcome_obj = outcome_obj,
                    eqtl_region = eqtl_region,
                    coloc_support = "Not_run",
                    note = collapse_notes(
                        if (!is.null(fetch_status)) fetch_status$note else NA_character_,
                        "No usable exact outcome rows before overlap step"
                    )
                )
                next
            }
            
            one <- run_single_target_coloc_abf(
                eqtl_region = eqtl_region,
                gwas_region = gwas_region,
                gene = g,
                outcome_obj = outcome_obj,
                protocol_cfg = protocol_cfg,
                run_cfg = run_cfg,
                preferred_df = preferred_df,
                fetch_status = fetch_status,
                err_file = err_file
            )
            
            res_list[[length(res_list) + 1L]] <- one
        }                                             
    }                                               
    
    dplyr::bind_rows(res_list)                     
}                                                