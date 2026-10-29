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
# The published file keeps every original reading, so core9/nutrient_chla5 are
# counted here before the analysis-time chlorophyll-a quality-control rule.
expected_counts <- c(514L, 391L, 440L, 459L, 439L, 454L)
chla_invalid <- sum(data$Chla_ugL < 0.01, na.rm = TRUE)
if (chla_invalid != 1L) stop("Unexpected number of sub-detection chlorophyll-a readings")
# After the rule is applied the analysis operates on 390 and 438 events.
analysis_data <- data
analysis_data$Chla_ugL[analysis_data$Chla_ugL < 0.01] <- NA
analysis_counts <- c(
  core9 = sum(complete.cases(analysis_data[, variables])),
  nutrient_chla5 = sum(complete.cases(analysis_data[, nutrient_chla5]))
)
if (!identical(as.integer(analysis_counts), c(390L, 438L))) stop("Analysis cohort counts changed")
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
  file.path("outputs", "results", "sessionInfo.txt"),
  # First-round revision outputs.
  file.path("outputs", "figures", "supplementary", "FigS6_som_map_layout_reproducibility.png"),
  file.path("outputs", "figures", "supplementary", "FigS7_transfer_matched_omission_all_k.png"),
  file.path("outputs", "tables", c(
    "Table4_decision_heuristic_r1.csv", "TableS3_som_grid_summary_r1_display.csv",
    "TableS11_survey_exclusion_r1_display.csv", "TableS15_sample_metrics_display.csv",
    "TableS16_jaccard_median_vs_mean_display.csv", "TableS17_som_map_reproducibility_display.csv",
    "TableS18_som_second_level_display.csv", "TableS19_ward_kmeans_allocation_display.csv",
    "TableS20_transfer_matched_omission_display.csv", "TableS21_modules_like_for_like_display.csv",
    "TableS_r1_resampling_training_membership_200.csv", "TableS_r1_resampling_seeds_200.csv",
    "TableS_r1_resample_som_refit_check.csv", "TableS_r1_missing_distance_identity_check.csv"
  )),
  file.path("outputs", "results", "results_facts_r1.json")
)
missing_outputs <- tracked_outputs[!file.exists(tracked_outputs)]
if (length(missing_outputs)) stop("Missing released outputs: ", paste(missing_outputs, collapse = ", "))

# The exported training memberships must describe 200 resamples of the 390
# primary events. Each resample keeps whole site-within-survey blocks and, in
# every survey, floor(0.8 x number of blocks) of them.
membership <- read.csv(file.path("outputs", "tables", "TableS_r1_resampling_training_membership_200.csv"),
                       check.names = FALSE, stringsAsFactors = FALSE)
resample_columns <- grep("^R[0-9]{3}$", names(membership), value = TRUE)
if (nrow(membership) != 390L || length(resample_columns) != 200L) stop("Unexpected resampling membership table")
if (!all(membership$event_id %in% data$event_id)) stop("Resampling membership refers to unknown events")
blocks <- unique(membership[, c("block_id_public", "campaign")])
expected_blocks <- floor(0.8 * table(blocks$campaign))
for (column in resample_columns) {
  within_block <- tapply(membership[[column]], membership$block_id_public, function(x) length(unique(x)))
  if (any(within_block != 1L)) stop("A resampling block was split in ", column)
  in_training <- tapply(membership[[column]], membership$block_id_public, max)
  kept <- tapply(in_training, blocks$campaign[match(names(in_training), blocks$block_id_public)], sum)
  if (!identical(as.integer(kept[names(expected_blocks)]), as.integer(expected_blocks))) {
    stop("Resample ", column, " does not keep 80% of the blocks in every survey")
  }
}

# Deterministic checks recorded by the revision scripts.
refit <- read.csv(file.path("outputs", "tables", "TableS_r1_resample_som_refit_check.csv"))
if (any(refit$mismatched_labels != 0)) stop("Refitted grouped-resample SOMs differ from the stored partitions")
identity <- read.csv(file.path("outputs", "tables", "TableS_r1_missing_distance_identity_check.csv"))
if (nrow(identity) != 50L || !all(identity$identical_partitions_k2_8) || any(identity$max_codebook_difference != 0)) {
  stop("Missing-distance SOM does not reduce to the 50 primary maps on complete events")
}
# Transfer counts use only held-out sets whose reference has more than one cluster.
transfer <- read.csv(file.path("outputs", "tables", "TableS_r1_transfer_vs_matched_omission.csv"))
needed <- c("loso_ami", "ami_lower", "reference_clusters", "reference_largest_share")
if (!all(needed %in% names(transfer)) || nrow(transfer) != 4L * 4L * 7L) stop("Unexpected transfer calibration table")
informative <- transfer$held_out_campaign != "2024-2025" & transfer$reference_clusters > 1L
facts <- jsonlite::fromJSON(file.path("outputs", "results", "results_facts_r1.json"))
tally <- facts$transfer_tallies$informative_surveys
if (tally$combinations != sum(informative) || tally$ari_below != sum(transfer$loso_ari[informative] < transfer$ari_lower[informative]) ||
    tally$ami_below != sum(transfer$loso_ami[informative] < transfer$ami_lower[informative])) {
  stop("Transfer tallies in results_facts_r1.json do not match the calibration table")
}

all_files <- list.files(".", recursive = TRUE, all.files = TRUE, full.names = TRUE)
all_files <- all_files[file.info(all_files)$isdir %in% FALSE]
# Prevent excluded cartographic and graphical-abstract assets from entering the release.
forbidden_names <- c("Fig1", "Fig2", "Graphical_Abstract", "study_map", "standard_map")
if (any(vapply(forbidden_names, function(x) any(grepl(x, all_files, fixed = TRUE)), logical(1)))) {
  stop("Excluded figure or map assets found in public release")
}

message("Release validation passed")
message("Complete-case counts in the released file (before the chlorophyll-a QC rule):")
print(counts)
message("Analysis cohorts after the QC rule (as reported in the manuscript):")
print(c(counts["all_events"], core9 = analysis_counts[["core9"]], counts[c("field4", "nutrient4")],
        nutrient_chla5 = analysis_counts[["nutrient_chla5"]], counts["partial5"]))
