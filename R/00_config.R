# Shared packages, paths, analysis settings and safe output helpers.
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(kohonen)
  library(mclust)
})

# Resolve the current script path for command-line and sourced execution.
script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit)) {
    candidate <- sub("^--file=", "", hit[1])
    if (nzchar(candidate) && candidate != "-" && file.exists(candidate)) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }
  frames <- sys.frames()
  source_files <- vapply(frames, function(frame) {
    value <- frame$ofile
    if (is.null(value)) NA_character_ else as.character(value)
  }, character(1))
  source_files <- source_files[!is.na(source_files) & file.exists(source_files)]
  if (length(source_files)) return(normalizePath(tail(source_files, 1), mustWork = TRUE))
  normalizePath(getwd(), mustWork = TRUE)
}

# Use repository-relative paths so the workflow is portable after cloning.
ANALYSIS_DIR <- dirname(script_path())
REPO_ROOT <- dirname(ANALYSIS_DIR)
DATA_FILE <- file.path(REPO_ROOT, "data", "analysis_events.csv")

RESULT_DIR <- file.path(REPO_ROOT, "outputs", "results")
MODEL_DIR <- file.path(RESULT_DIR, "models")
TABLE_DIR <- file.path(REPO_ROOT, "outputs", "tables")
FIG_MAIN_DIR <- file.path(REPO_ROOT, "outputs", "figures", "main")
FIG_SUPP_DIR <- file.path(REPO_ROOT, "outputs", "figures", "supplementary")

invisible(lapply(
  c(RESULT_DIR, MODEL_DIR, TABLE_DIR, FIG_MAIN_DIR, FIG_SUPP_DIR),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

# Define the common feature space and prespecified variable modules.
VARS <- c(
  "pH", "WT_C", "DO_mgL", "Cond_uScm", "Chla_ugL",
  "TN_mgL", "NO3_mgL", "TP_mgL", "SRP_mgL"
)
LOG_VARS <- c("Chla_ugL", "TN_mgL", "NO3_mgL", "TP_mgL", "SRP_mgL")
FIELD_VARS <- c("pH", "WT_C", "DO_mgL", "Cond_uScm")
NUTRIENT_VARS <- c("TN_mgL", "NO3_mgL", "TP_mgL", "SRP_mgL")
NUTRIENT_CHLA_VARS <- c("Chla_ugL", NUTRIENT_VARS)
CAMPAIGNS <- c("2012-2013", "2017-2018", "2020-2021", "2024-2025")
K_VALUES <- 2:8
GRIDS <- data.frame(
  xdim = c(5L, 7L, 8L, 10L, 12L),
  ydim = c(5L, 5L, 6L, 7L, 8L)
)
# Fix the model grids, random seeds and resampling budget used in the manuscript.
SOM_SEEDS <- 1101:1150
SOM_RLEN <- 1000L
N_RESAMPLE <- 200L
RESAMPLE_FRACTION <- 0.80
RESAMPLE_SEEDS <- 71001 + seq_len(N_RESAMPLE)
TRANSFER_SEEDS <- 81001:81010

# Quality-control threshold applied in 02_prepare_data.R. Chlorophyll-a readings
# below this value are two orders of magnitude under any routine determination
# limit and are treated as invalid rather than as extremely low concentrations.
# The published data file keeps the original reading; only the analysis excludes it.
CHLA_MIN_VALID <- 0.01

# Fixed grid used in 04b_fixed_grid_som.R for the like-for-like stability
# comparison against the baselines, which have no grid hyperparameter. 8x6 has
# the highest ten-neighbour preservation in the grid sweep. The initialisation
# seed is held constant so that the remaining variation comes from sample
# composition alone, as it does for the deterministic baselines.
FIXED_GRID <- c(xdim = 8L, ydim = 6L)
FIXED_GRID_SEED <- 91001L

# Outlier sensitivity in 06b_outlier_sensitivity.R: events whose standardised
# value on any variable exceeds this magnitude are dropped to test whether heavy
# tails drive the fragile-cluster result.
OUTLIER_Z_THRESHOLD <- 4

# Guard against accidental changes in the released analytical cohorts.
EXPECTED_COUNTS <- c(
  all_surface = 514L,
  core9 = 390L,
  field4 = 440L,
  nutrient4 = 459L,
  nutrient_chla5 = 438L,
  partial5 = 454L
)

# Create parent directories before writing reproducible outputs.
# First-round revision checks (scripts 10-14, October 2026). Every SOM
# sensitivity ensemble uses the five grids crossed with the first ten SOM seeds,
# as in 06_sensitivity.R. Matched-size omission replicates remove randomly chosen
# site-within-survey blocks from all surveys until the omitted event count is
# closest to the size of one held-out survey.
SENSITIVITY_SEEDS <- SOM_SEEDS[1:10]
MATCHED_OMISSION_REPLICATES <- 50L
MATCHED_OMISSION_SEEDS <- 95000L + seq_len(MATCHED_OMISSION_REPLICATES)
EXCLUSION_OMISSION_SEEDS <- 96000L + seq_len(MATCHED_OMISSION_REPLICATES)

write_csv_safe <- function(x, path) {
  readr::write_excel_csv(x, path, na = "")
  invisible(path)
}

save_rds_safe <- function(x, path) {
  saveRDS(x, path, version = 3, compress = "xz")
  invisible(path)
}

# Load functions shared by all analysis stages.
source(file.path(ANALYSIS_DIR, "01_helpers.R"), local = FALSE)
