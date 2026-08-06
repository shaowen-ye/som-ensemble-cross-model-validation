args <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
repo_root <- dirname(normalizePath(this_file, mustWork = TRUE))
analysis_dir <- file.path(repo_root, "R")
rscript <- Sys.which("Rscript")
if (!nzchar(rscript)) stop("Rscript is not available on PATH")

pipeline <- c(
  "02_prepare_data.R",
  "03_fit_som_ensemble.R",
  "04_resampling_validation.R",
  "05_campaign_transfer.R",
  "06_sensitivity.R",
  "07_gradient_robustness.R",
  "08_summarise_results.R",
  "09_make_figures.R"
)

for (script in pipeline) {
  path <- file.path(analysis_dir, script)
  message("\n==> ", script)
  status <- system2(rscript, path)
  if (!identical(status, 0L)) stop("Pipeline stopped at ", script, " (status ", status, ")")
}

message("\nAnalysis rebuild complete. See outputs/ for figures, tables, and session information.")
