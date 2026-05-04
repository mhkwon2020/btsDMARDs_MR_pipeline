# R/00_utils.R
# Core utility helpers
# Version: 1.0
# Author: Mihye Kwon
suppressPackageStartupMessages({
  library(fs)
  library(glue)
  library(yaml)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

log_message <- function(log_file, ...) {
  cat(
    sprintf(
      "[%s] %s\n",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      glue(..., .envir = parent.frame())
    ),
    file = log_file,
    append = TRUE
  )
}

write_yaml_snapshot <- function(obj, path) {
  yaml::write_yaml(obj, path)
}

write_csv_safe <- function(df, path) {
  if (is.null(df)) df <- data.frame()
  utils::write.csv(df, path, row.names = FALSE, na = "")
}

make_run_dirs <- function(base_dir = "runs", prefix = "eqtl_abstract") {
  ts <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  run_dir <- fs::path(base_dir, paste0(prefix, "_", ts))
  results_dir <- fs::path(run_dir, "results")
  figures_dir <- fs::path(results_dir, "figures")
  artifacts_dir <- fs::path(run_dir, "artifacts")

  fs::dir_create(results_dir, recurse = TRUE)
  fs::dir_create(figures_dir, recurse = TRUE)
  fs::dir_create(artifacts_dir, recurse = TRUE)

  list(
    run_dir = run_dir,
    results_dir = results_dir,
    figures_dir = figures_dir,
    artifacts_dir = artifacts_dir,
    log_file = fs::path(run_dir, "pipeline_log.txt"),
    err_file = fs::path(run_dir, "pipeline_errors.txt")
  )
}

save_session_info <- function(path) {
  sink(path)
  print(sessionInfo())
  sink()
}

preferred_method_order <- function(run_cfg) {
  run_cfg$mr_methods$preferred_method_order %||%
    c("Inverse variance weighted", "Wald ratio", "Weighted median", "MR Egger")
}

pick_preferred_result <- function(df, method_order) {
  if (is.null(df) || nrow(df) == 0) return(df[0, , drop = FALSE])
  df$.method_rank <- match(df$method, method_order)
  df$.method_rank[is.na(df$.method_rank)] <- 999L
  df <- df[order(df$.method_rank, df$pval), , drop = FALSE]
  df[1, setdiff(names(df), ".method_rank"), drop = FALSE]
}
