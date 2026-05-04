# R/14_additional_figures.R
# Additional triangulation figure utilities
# Version: 1.0
# Author: Mihye Kwon
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(fs)
})

fig_pick_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

fig_upper_gene <- function(x) {
  toupper(as.character(x))
}

fig_pretty_status <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- "Not available"

  out <- dplyr::case_when(
    x == "Strong_coloc" ~ "Strong coloc",
    x == "Distinct_signals_likely" ~ "Distinct signals",
    x == "Weak_or_uncertain" ~ "Weak/uncertain",
    x == "Strong" ~ "Strong",
    x == "Weak_or_partial" ~ "Weak/partial",
    x == "Not_supported" ~ "Not supported",
    x == "Non_estimable" ~ "Not estimable",
    x == "Supportive_primary_and_replication" ~ "Supportive",
    x == "Supportive_primary_only" ~ "Primary only",
    x == "No_support_detected" ~ "No support",
    x == "concordant_support_two_source" ~ "Concordant",
    x == "discordant_support_two_source" ~ "Discordant",
    x == "two_source_estimable_no_dual_support" ~ "2-source estimable",
    x == "single_source_only" ~ "Single source",
    x == "source_unavailable_both" ~ "Unavailable",
    x == "not_assessable" ~ "Not assessable",
    x == "NA" ~ "Not available",
    TRUE ~ gsub("_", " ", x, fixed = TRUE)
  )
  out
}

fig_status_group <- function(x) {
  dplyr::case_when(
    x %in% c("Strong_coloc", "Strong", "Supportive_primary_and_replication", "concordant_support_two_source") ~ "Strong support",
    x %in% c("Weak_or_partial", "Supportive_primary_only", "Distinct_signals_likely", "Weak_or_uncertain", "two_source_estimable_no_dual_support", "single_source_only") ~ "Intermediate",
    x %in% c("Not_supported", "No_support_detected", "discordant_support_two_source") ~ "No support",
    x %in% c("Non_estimable", "Not_estimable", "source_unavailable_both", "not_assessable", "NA") ~ "Missing / not estimable",
    TRUE ~ "Other"
  )
}

tri_build_heatmap_df_v2 <- function(master_table) {
  class_col <- fig_pick_col(master_table, c("class.y", "class", "benchmark_class"))
  coloc_col <- fig_pick_col(master_table, c("final_coloc_support", "coloc_support.y", "coloc_support.x", "coloc_support"))

  out <- master_table %>%
    dplyr::transmute(
      gene = fig_upper_gene(.data$gene),
      `eQTL evidence` = .data[[class_col]],
      `Colocalisation` = .data[[coloc_col]],
      `UKB-PPP pQTL` = .data$pqtl_ukbppp_support %||% NA_character_,
      `deCODE pQTL` = .data$pqtl_decode_support %||% NA_character_,
      `Cross-source robustness` = .data$cross_source_robustness %||% NA_character_
    ) %>%
    tidyr::pivot_longer(cols = -gene, names_to = "layer", values_to = "status") %>%
    dplyr::mutate(
      layer = factor(layer, levels = c(
        "eQTL evidence",
        "Colocalisation",
        "UKB-PPP pQTL",
        "deCODE pQTL",
        "Cross-source robustness"
      )),
      status_raw = ifelse(is.na(status) | !nzchar(status), "NA", status),
      status_label = fig_pretty_status(status_raw),
      status_group = fig_status_group(status_raw)
    )
  out
}

write_additional_triangulation_figures_v2 <- function(dirs,
                                                      master_table,
                                                      pqtl_source_summary = NULL,
                                                      priority_coloc_df = NULL,
                                                      overlap_dt = NULL,
                                                      log_file = NULL) {
  if (is.null(master_table) || !is.data.frame(master_table) || nrow(master_table) == 0) {
    return(invisible(NULL))
  }

  fs::dir_create(dirs$figures_dir, recurse = TRUE)
  gene_order <- fig_upper_gene(master_table$gene)

  # Fig 6
  heat_df <- tri_build_heatmap_df_v2(master_table) %>%
    dplyr::mutate(gene = factor(gene, levels = rev(gene_order)))

  p_heat <- ggplot(heat_df, aes(x = layer, y = gene, fill = status_group)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = status_label), size = 3.0, lineheight = 0.9) +
    scale_fill_manual(
      values = c(
        "Strong support" = "#1b9e77",
        "Intermediate" = "#d95f02",
        "No support" = "#7570b3",
        "Missing / not estimable" = "#bdbdbd",
        "Other" = "#e6ab02"
      ),
      drop = FALSE,
      name = "Evidence group"
    ) +
    labs(
      title = "Triangulation summary across evidence layers",
      subtitle = "Evidence progresses from eQTL and locus support to source-specific pQTL robustness",
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1),
      axis.text.y = element_text(face = "bold"),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 11),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(10, 15, 10, 10)
    )

  ggsave(fs::path(dirs$figures_dir, "Fig6_triangulation_heatmap.png"), p_heat,
         width = 12.8, height = 8.8, dpi = 300, bg = "white")

  # Fig 7
  src_df <- master_table %>%
    dplyr::select(gene, pqtl_ukbppp_support, pqtl_decode_support) %>%
    tidyr::pivot_longer(
      cols = c(pqtl_ukbppp_support, pqtl_decode_support),
      names_to = "source",
      values_to = "support"
    ) %>%
    dplyr::mutate(
      gene = factor(fig_upper_gene(gene), levels = rev(gene_order)),
      source = dplyr::recode(source,
                             pqtl_ukbppp_support = "UKB-PPP",
                             pqtl_decode_support = "deCODE"),
      source = factor(source, levels = c("UKB-PPP", "deCODE")),
      support_raw = ifelse(is.na(support) | !nzchar(support), "Not_estimable", support),
      support_label = dplyr::case_when(
        support_raw == "Supportive_primary_and_replication" ~ "Supportive",
        support_raw == "Supportive_primary_only" ~ "Primary only",
        support_raw == "No_support_detected" ~ "No support",
        support_raw == "Not_estimable" ~ "Not estimable",
        TRUE ~ gsub("_", " ", support_raw, fixed = TRUE)
      ),
      support_group = dplyr::case_when(
        support_raw == "Supportive_primary_and_replication" ~ "Supportive",
        support_raw == "Supportive_primary_only" ~ "Partial support",
        support_raw == "No_support_detected" ~ "No support",
        TRUE ~ "Not estimable"
      )
    )

  p_src <- ggplot(src_df, aes(x = source, y = gene, fill = support_group)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = support_label), size = 3.2) +
    scale_fill_manual(
      values = c(
        "Supportive" = "#4C84E0",
        "Partial support" = "#7DB7FF",
        "No support" = "#F1786B",
        "Not estimable" = "#00BA38"
      ),
      drop = FALSE,
      name = "pQTL support"
    ) +
    labs(
      title = "Source-specific pQTL support",
      subtitle = "UKB-PPP and deCODE compared at the gene level",
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      axis.text.y = element_text(face = "bold"),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 11),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    )

  ggsave(fs::path(dirs$figures_dir, "Fig7_pqtl_source_support_matrix.png"), p_src,
         width = 7.4, height = 9.0, dpi = 300, bg = "white")

  # Fig 8
  if (!is.null(priority_coloc_df) && is.data.frame(priority_coloc_df) && nrow(priority_coloc_df) > 0) {
    coloc_support_col <- fig_pick_col(priority_coloc_df, c("final_coloc_support", "coloc_support"))
    pc <- priority_coloc_df %>%
      dplyr::filter(outcome_name == "RA_EUR_PRIMARY") %>%
      dplyr::mutate(
        gene = factor(fig_upper_gene(gene), levels = fig_upper_gene(gene)),
        support_label = fig_pretty_status(.data[[coloc_support_col]])
      )

    p_coloc <- ggplot(pc, aes(x = gene, y = final_PP.H4)) +
      geom_segment(aes(xend = gene, y = 0, yend = final_PP.H4), linewidth = 0.6) +
      geom_point(size = 4) +
      geom_text(aes(label = support_label), nudge_y = 0.05, size = 4, check_overlap = TRUE) +
      coord_cartesian(ylim = c(0, 1.08)) +
      labs(
        title = "Priority-target colocalisation summary",
        x = NULL,
        y = "Final PP.H4"
      ) +
      theme_bw(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold")
      )

    ggsave(fs::path(dirs$figures_dir, "Fig8_priority_colocalisation_summary.png"), p_coloc,
           width = 8.5, height = 5.2, dpi = 300, bg = "white")
  }

  # Fig S1
  if (!is.null(overlap_dt) && is.data.frame(overlap_dt) && nrow(overlap_dt) > 0) {
    ov <- overlap_dt %>%
      dplyr::mutate(
        exposure_label = dplyr::recode(
          exposure_id,
          EQTLGEN_BLOOD = "eQTLGen blood cis-eQTL",
          UKB_PPP_EUR = "UKB-PPP (EUR)",
          deCODE_EUR = "deCODE (EUR)",
          .default = exposure_name
        ),
        outcome_label = dplyr::recode(
          outcome_name,
          RA_EUR_PRIMARY = "RA EUR primary",
          RA_EUR_REP1_FinnGen = "RA EUR replication",
          RA_EUR_SEROPOS = "RA seropositive",
          RA_EUR_SERONEG = "RA seronegative",
          RA_EAS_BBJ = "RA EAS BBJ",
          .default = outcome_name
        ),
        exposure_label = factor(exposure_label, levels = rev(unique(exposure_label))),
        outcome_label = factor(outcome_label, levels = c(
          "RA EUR primary", "RA EUR replication", "RA seropositive", "RA seronegative", "RA EAS BBJ"
        )),
        overlap_group = dplyr::case_when(
          overlap_status == "possible_but_unquantified_overlap" ~ "Possible overlap",
          overlap_status == "no_detected_direct_cohort_overlap" ~ "No direct overlap",
          TRUE ~ "Other / unknown"
        ),
        overlap_label = dplyr::case_when(
          overlap_group == "Possible overlap" ~ "Possible",
          overlap_group == "No direct overlap" ~ "No direct",
          TRUE ~ "Other"
        )
      )

    p_overlap <- ggplot(ov, aes(x = outcome_label, y = exposure_label, fill = overlap_group)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = overlap_label), size = 3.3) +
      scale_fill_manual(
        values = c(
          "No direct overlap" = "#F1786B",
          "Possible overlap" = "#12B7BD",
          "Other / unknown" = "#BDBDBD"
        ),
        drop = FALSE,
        name = "Overlap status"
      ) +
      labs(
        title = "Sample overlap screening matrix",
        x = NULL,
        y = NULL
      ) +
      theme_bw(base_size = 13) +
      theme(
        axis.text.x = element_text(angle = 32, hjust = 1, vjust = 1),
        axis.text.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")
      )

    ggsave(fs::path(dirs$figures_dir, "FigS1_sample_overlap_matrix.png"), p_overlap,
           width = 8.8, height = 4.3, dpi = 300, bg = "white")
  }

  invisible(NULL)
}


# Backward-compatible wrapper
write_additional_triangulation_figures <- function(dirs,
                                                   master_table,
                                                   pqtl_source_summary = NULL,
                                                   priority_coloc_df = NULL,
                                                   overlap_dt = NULL,
                                                   log_file = NULL) {
  write_additional_triangulation_figures_v2(
    dirs = dirs,
    master_table = master_table,
    pqtl_source_summary = pqtl_source_summary,
    priority_coloc_df = priority_coloc_df,
    overlap_dt = overlap_dt,
    log_file = log_file
  )
}
