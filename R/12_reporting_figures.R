# R/12_reporting_figures.R
# Reporting figure utilities
# Version: 1.0
# Author: Mihye Kwon
suppressPackageStartupMessages({
  library(dplyr)
})

write_full_study_figures <- function(dirs,
                                     protocol_cfg,
                                     attrition_df,
                                     preferred_df,
                                     class_support_df,
                                     cross_proxy_df = NULL,
                                     log_file = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    if (!is.null(log_file)) log_message(log_file, "ggplot2 not available; figures skipped")
    return(invisible(NULL))
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    if (!is.null(log_file)) log_message(log_file, "tidyr not available; figures skipped")
    return(invisible(NULL))
  }

  gene_levels <- tolower(protocol_cfg$targets$gene_list)

  attr_plot <- attrition_df %>%
    group_by(gene) %>%
    summarise(
      raw = suppressWarnings(max(raw_variant_count, na.rm = TRUE)),
      post_p = suppressWarnings(max(post_p_threshold_count, na.rm = TRUE)),
      post_clump = suppressWarnings(max(post_clump_count, na.rm = TRUE)),
      usable = suppressWarnings(max(post_harmonisation_usable_count, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(
      cols = c(raw, post_p, post_clump, usable),
      names_to = "stage",
      values_to = "n"
    ) %>%
    mutate(gene = factor(gene, levels = gene_levels))

  p1 <- ggplot2::ggplot(attr_plot, ggplot2::aes(x = gene, y = n, fill = stage)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::labs(title = "Instrument attrition overview", x = NULL, y = "Count") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  ggplot2::ggsave(
    fs::path(dirs$figures_dir, "Fig1_attrition_overview.png"),
    p1, width = 12, height = 6, dpi = 300
  )

  prim <- preferred_df %>%
    filter(outcome_name == "RA_EUR_PRIMARY") %>%
    select(gene, b_primary = b, p_primary = pval)

  rep1 <- preferred_df %>%
    filter(outcome_name == "RA_EUR_REP1_FinnGen") %>%
    select(gene, b_rep = b, p_rep = pval)

  sc <- full_join(prim, rep1, by = "gene")

  if (nrow(sc) > 0) {
    p2 <- ggplot2::ggplot(sc, ggplot2::aes(x = b_primary, y = b_rep, label = gene)) +
      ggplot2::geom_point() +
      ggplot2::geom_hline(yintercept = 0, linetype = 2) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2) +
      ggplot2::geom_text(nudge_y = 0.01, check_overlap = TRUE) +
      ggplot2::labs(
        title = "Primary vs replication effect estimates",
        x = "Primary beta",
        y = "Replication beta"
      ) +
      ggplot2::theme_bw()

    ggplot2::ggsave(
      fs::path(dirs$figures_dir, "Fig2_primary_replication_scatter.png"),
      p2, width = 7, height = 6, dpi = 300
    )
  }

  if (nrow(class_support_df) > 0) {
    heat <- class_support_df %>%
      mutate(
        gene = factor(gene, levels = rev(gene_levels)),
        class = factor(.data$class, levels = c("Strong", "Weak_or_partial", "Not_supported", "Non_estimable"))
      )

    p3 <- ggplot2::ggplot(heat, ggplot2::aes(x = "RA benchmark", y = gene, fill = class)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::labs(title = "Benchmark support summary", x = NULL, y = NULL) +
      ggplot2::theme_bw()

    ggplot2::ggsave(
      fs::path(dirs$figures_dir, "Fig3_benchmark_heatmap.png"),
      p3, width = 5, height = 8, dpi = 300
    )
  }

  forest_df <- preferred_df %>%
    filter(outcome_name == "RA_EUR_PRIMARY") %>%
    arrange(pval)

  if (nrow(forest_df) > 0) {
    forest_df$gene <- factor(forest_df$gene, levels = rev(gene_levels[gene_levels %in% forest_df$gene]))

    p4 <- ggplot2::ggplot(
      forest_df,
      ggplot2::aes(x = gene, y = b, ymin = b - 1.96 * se, ymax = b + 1.96 * se)
    ) +
      ggplot2::geom_pointrange() +
      ggplot2::geom_hline(yintercept = 0, linetype = 2) +
      ggplot2::coord_flip() +
      ggplot2::labs(title = "Primary benchmark estimates", x = NULL, y = "Beta (95% CI)") +
      ggplot2::theme_bw()

    ggplot2::ggsave(
      fs::path(dirs$figures_dir, "Fig4_selected_forest.png"),
      p4, width = 8, height = 7, dpi = 300
    )
  }

  if (!is.null(cross_proxy_df) && nrow(cross_proxy_df) > 0) {
    proxy_plot <- cross_proxy_df %>%
      mutate(
        gene = factor(gene, levels = rev(gene_levels)),
        cross_proxy_consistency = factor(
          cross_proxy_consistency,
          levels = c("Concordant", "Mixed", "Discordant", "No_pQTL_data")
        )
      )

    p5 <- ggplot2::ggplot(proxy_plot, ggplot2::aes(x = "eQTL vs pQTL", y = gene, fill = cross_proxy_consistency)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::labs(title = "Cross-proxy consistency", x = NULL, y = NULL) +
      ggplot2::theme_bw()

    ggplot2::ggsave(
      fs::path(dirs$figures_dir, "Fig5_cross_proxy_consistency.png"),
      p5, width = 5.5, height = 8, dpi = 300
    )
  }

  invisible(NULL)
}
