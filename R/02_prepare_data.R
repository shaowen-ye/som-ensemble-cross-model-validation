args <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "00_config.R"))

assert_packages(c("jsonlite", "kohonen", "mclust", "cluster"))
if (!file.exists(DATA_FILE)) stop("Analysis-ready input is unavailable: ", DATA_FILE)

raw <- readr::read_csv(DATA_FILE, show_col_types = FALSE)
required <- c(
  "event_id", "campaign", "block_id", VARS
)
missing_columns <- setdiff(required, names(raw))
if (length(missing_columns)) stop("Missing required fields: ", paste(missing_columns, collapse = ", "))

data <- raw %>%
  transmute(
    record_id = event_id,
    sampling_event_id = event_id,
    campaign = factor(campaign, levels = CAMPAIGNS),
    block_id,
    across(all_of(VARS))
  )

if (anyDuplicated(data$record_id) || anyDuplicated(data$sampling_event_id)) {
  stop("Surface event identifiers are not unique")
}
if (any(is.na(data$campaign))) stop("Unexpected campaign label")

counts <- c(
  all_surface = nrow(data),
  core9 = sum(complete.cases(data[, VARS])),
  field4 = sum(complete.cases(data[, FIELD_VARS])),
  nutrient4 = sum(complete.cases(data[, NUTRIENT_VARS])),
  nutrient_chla5 = sum(complete.cases(data[, NUTRIENT_CHLA_VARS])),
  partial5 = sum(rowSums(!is.na(data[, VARS])) >= 5)
)
if (!identical(as.integer(counts[names(EXPECTED_COUNTS)]), as.integer(EXPECTED_COUNTS))) {
  stop("Frozen cohort counts changed: ", paste(names(counts), counts, collapse = "; "))
}

cohorts <- list(
  all_surface = data,
  core9 = data[complete.cases(data[, VARS]), ],
  field4 = data[complete.cases(data[, FIELD_VARS]), ],
  nutrient4 = data[complete.cases(data[, NUTRIENT_VARS]), ],
  nutrient_chla5 = data[complete.cases(data[, NUTRIENT_CHLA_VARS]), ],
  partial5 = data[rowSums(!is.na(data[, VARS])) >= 5, ]
)

fit <- transform_fit(cohorts$core9)
x <- transform_apply(cohorts$core9, fit)
if (any(!is.finite(x))) stop("Primary transformed matrix contains non-finite values")

manifest <- list(
  input_file = "data/analysis_events.csv",
  input_rows = nrow(raw),
  input_columns = ncol(raw),
  counts = as.list(counts),
  variables = VARS,
  log10_variables = LOG_VARS,
  primary_independence_unit = "unique surface sampling event",
  resampling_block = "site within campaign"
)

save_rds_safe(list(data = data, cohorts = cohorts, transform = fit, x = x, manifest = manifest), file.path(RESULT_DIR, "prepared_data.rds"))
jsonlite::write_json(manifest, file.path(RESULT_DIR, "input_manifest.json"), auto_unbox = TRUE, pretty = TRUE)

count_table <- bind_rows(lapply(names(cohorts), function(nm) {
  cohorts[[nm]] %>% dplyr::count(campaign, name = "n") %>% mutate(cohort = nm, .before = 1)
}))
write_csv_safe(count_table, file.path(TABLE_DIR, "TableS_cohort_counts.csv"))

variable_table <- bind_rows(lapply(VARS, function(v) {
  z <- data[[v]]
  data.frame(
    variable = v,
    available_n = sum(!is.na(z)),
    missing_n = sum(is.na(z)),
    minimum = min(z, na.rm = TRUE),
    median = median(z, na.rm = TRUE),
    maximum = max(z, na.rm = TRUE),
    primary_transformation = ifelse(v %in% LOG_VARS, "log10 then z-score", "z-score")
  )
}))
write_csv_safe(variable_table, file.path(TABLE_DIR, "Table1_variables_and_transformations.csv"))

variable_labels <- c(
  pH = "pH", WT_C = "Water temperature", DO_mgL = "Dissolved oxygen",
  Cond_uScm = "Conductivity", Chla_ugL = "Chlorophyll-a", TN_mgL = "Total nitrogen",
  NO3_mgL = "Nitrate nitrogen", TP_mgL = "Total phosphorus", SRP_mgL = "Soluble reactive phosphorus"
)
variable_units <- c(
  pH = "unitless", WT_C = "°C", DO_mgL = "mg L−1", Cond_uScm = "µS cm−1",
  Chla_ugL = "µg L−1", TN_mgL = "mg L−1", NO3_mgL = "mg L−1",
  TP_mgL = "mg L−1", SRP_mgL = "mg L−1"
)
variable_table_display <- variable_table %>%
  transmute(
    Variable = unname(variable_labels[variable]), Unit = unname(variable_units[variable]),
    `Available n` = available_n, `Missing n` = missing_n,
    Minimum = signif(minimum, 3), Median = signif(median, 3), Maximum = signif(maximum, 3),
    `Primary transformation` = ifelse(variable %in% LOG_VARS, "log10, then z standardisation", "z standardisation")
  )
write_csv_safe(variable_table_display, file.path(TABLE_DIR, "Table1_variables_and_transformations_display.csv"))

cat("Prepared authoritative cohorts\n")
print(counts)
