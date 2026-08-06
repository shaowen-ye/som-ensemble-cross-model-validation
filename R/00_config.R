suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(kohonen)
  library(mclust)
})

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
SOM_SEEDS <- 1101:1150
SOM_RLEN <- 1000L
N_RESAMPLE <- 200L
RESAMPLE_FRACTION <- 0.80
RESAMPLE_SEEDS <- 71001 + seq_len(N_RESAMPLE)
TRANSFER_SEEDS <- 81001:81010

EXPECTED_COUNTS <- c(
  all_surface = 514L,
  core9 = 391L,
  field4 = 440L,
  nutrient4 = 459L,
  nutrient_chla5 = 439L,
  partial5 = 454L
)

write_csv_safe <- function(x, path) {
  readr::write_excel_csv(x, path, na = "")
  invisible(path)
}

save_rds_safe <- function(x, path) {
  saveRDS(x, path, version = 3, compress = "xz")
  invisible(path)
}

source(file.path(ANALYSIS_DIR, "01_helpers.R"), local = FALSE)
