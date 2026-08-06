# Validate the public data boundary, analytical cohort counts and released outputs.
data_file <- file.path("data", "analysis_events.csv")
if (!file.exists(data_file)) stop("Missing public analysis data")

data <- read.csv(data_file, check.names = FALSE, na.strings = "", stringsAsFactors = FALSE)
variables <- c(
  "pH", "WT_C", "DO_mgL", "Cond_uScm", "Chla_ugL",
  "TN_mgL", "NO3_mgL", "TP_mgL", "SRP_mgL"
)
expected_columns <- c("event_id", "campaign", "block_id", variables)
if (!identical(names(data), expected_columns)) stop("Unexpected public data columns")
if (nrow(data) != 514L || ncol(data) != 12L) stop("Unexpected public data dimensions")
if (anyDuplicated(data$event_id)) stop("Duplicated event identifiers")
if (anyDuplicated(data[, c("event_id", "campaign")])) stop("Duplicated events")
if (!identical(sort(unique(data$campaign)), sort(c("2012-2013", "2017-2018", "2020-2021", "2024-2025")))) {
  stop("Unexpected campaign labels")
}

field4 <- c("pH", "WT_C", "DO_mgL", "Cond_uScm")
nutrient4 <- c("TN_mgL", "NO3_mgL", "TP_mgL", "SRP_mgL")
nutrient_chla5 <- c("Chla_ugL", nutrient4)
# Recalculate the six sample-size statements used in the manuscript.
counts <- c(
  all_events = nrow(data),
  core9 = sum(complete.cases(data[, variables])),
  field4 = sum(complete.cases(data[, field4])),
  nutrient4 = sum(complete.cases(data[, nutrient4])),
  nutrient_chla5 = sum(complete.cases(data[, nutrient_chla5])),
  partial5 = sum(rowSums(!is.na(data[, variables])) >= 5)
)
expected_counts <- c(514L, 391L, 440L, 459L, 439L, 454L)
if (!identical(as.integer(counts), expected_counts)) stop("Cohort counts changed")

log_variables <- c("Chla_ugL", "TN_mgL", "NO3_mgL", "TP_mgL", "SRP_mgL")
if (any(vapply(data[log_variables], function(x) any(x <= 0, na.rm = TRUE), logical(1)))) {
  stop("Non-positive observations found in a log10 variable")
}

# Confirm that the repository contains the selected core result tables.
tracked_outputs <- c(
  file.path("outputs", "figures", "main", paste0("Fig", 3:6, c(
    "_som_components_and_topology.png",
    "_resampling_and_cross_model_stability.png", "_continuous_geometry_and_uncertainty.png",
    "_transfer_and_sensitivity.png"
  ))),
  file.path("outputs", "tables", "Table1_variables_and_transformations_display.csv"),
  file.path("outputs", "tables", "Table2_controlled_model_comparison.csv"),
  file.path("outputs", "tables", "Table3_resampling_stability_k2_k3.csv"),
  file.path("outputs", "tables", "Table4A_campaign_transfer_k3.csv"),
  file.path("outputs", "tables", "Table4B_sensitivity_k3.csv"),
  file.path("outputs", "tables", "Table4C_source_uncertainty_exclusion_k3.csv"),
  file.path("outputs", "results", "sessionInfo.txt")
)
missing_outputs <- tracked_outputs[!file.exists(tracked_outputs)]
if (length(missing_outputs)) stop("Missing released outputs: ", paste(missing_outputs, collapse = ", "))

all_files <- list.files(".", recursive = TRUE, all.files = TRUE, full.names = TRUE)
all_files <- all_files[file.info(all_files)$isdir %in% FALSE]
# Prevent excluded cartographic and graphical-abstract assets from entering the release.
forbidden_names <- c("Fig1", "Fig2", "Graphical_Abstract", "study_map", "standard_map")
if (any(vapply(forbidden_names, function(x) any(grepl(x, all_files, fixed = TRUE)), logical(1)))) {
  stop("Excluded figure or map assets found in public release")
}

message("Release validation passed")
print(counts)
