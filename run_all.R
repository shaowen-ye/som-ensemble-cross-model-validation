# Run the complete public workflow in manuscript order.
args <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
repo_root <- dirname(normalizePath(this_file, mustWork = TRUE))
analysis_dir <- file.path(repo_root, "R")
rscript <- Sys.which("Rscript")
if (!nzchar(rscript)) stop("Rscript is not available on PATH")

# Figures 1-2 and the graphical abstract are intentionally outside this release.
pipeline <- c(
  "02_prepare_data.R",
  "03_fit_som_ensemble.R",
  "04_resampling_validation.R",
  "04b_fixed_grid_som.R",
  "05_campaign_transfer.R",
  "06_sensitivity.R",
  "06b_outlier_sensitivity.R",
  "07_gradient_robustness.R",
  "08_summarise_results.R",
  # First-round revision checks; the figures below read their outputs.
  "10_r1_som_topology.R",
  "11_r1_partition_stages.R",
  "12_r1_metric_checks.R",
  "13_r1_transfer_calibration.R",
  "14_r1_sensitivity_like_for_like.R",
  "09_make_figures.R",
  "15_r1_make_tables.R"
)

# Stop immediately if any analytical stage fails.
for (script in pipeline) {
  path <- file.path(analysis_dir, script)
  message("\n==> ", script)
  status <- system2(rscript, path)
  if (!identical(status, 0L)) stop("Pipeline stopped at ", script, " (status ", status, ")")
}

message("\nAnalysis rebuild complete. See outputs/ for figures, tables, and session information.")
