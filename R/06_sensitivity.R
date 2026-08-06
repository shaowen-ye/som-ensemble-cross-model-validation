args <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "00_config.R"))

assert_packages(c("kohonen", "mclust"))
prepared <- readRDS(file.path(RESULT_DIR, "prepared_data.rds"))
resampling_file <- Sys.getenv("EI_RESAMPLING_FILE", unset = file.path(MODEL_DIR, "cross_model_resampling_200.rds"))
if (!file.exists(resampling_file)) stop("Resampling result is unavailable")
resampling <- readRDS(resampling_file)

fit_som_sensitivity_ensemble <- function(data, vars, log_vars, name, allow_partial = FALSE) {
  transform <- transform_fit(data, vars = vars, log_vars = log_vars)
  x <- transform_apply(data, transform)
  seeds <- SOM_SEEDS[1:10]
  n_models <- nrow(GRIDS) * length(seeds)
  assignments <- setNames(lapply(K_VALUES, function(k) matrix(NA_integer_, nrow(data), n_models)), paste0("k", K_VALUES))
  metric_rows <- list(); counter <- 0L
  complete_rows <- which(complete.cases(x))
  for (g in seq_len(nrow(GRIDS))) for (seed in seeds) {
    counter <- counter + 1L
    set.seed(seed)
    init <- NULL
    if (allow_partial) {
      nodes <- GRIDS$xdim[g] * GRIDS$ydim[g]
      init <- x[sample(complete_rows, nodes, replace = FALSE), , drop = FALSE]
    }
    fit_args <- list(
      X = x,
      grid = kohonen::somgrid(GRIDS$xdim[g], GRIDS$ydim[g], topo = "hexagonal"),
      rlen = SOM_RLEN,
      alpha = c(0.05, 0.01),
      keep.data = TRUE,
      mode = "online"
    )
    if (allow_partial) {
      fit_args$maxNA.fraction <- 4 / 9
      fit_args$init <- init
    }
    model <- do.call(kohonen::som, fit_args)
    metric_rows[[counter]] <- data.frame(
      analysis = name,
      model = counter,
      grid = sprintf("%dx%d", GRIDS$xdim[g], GRIDS$ydim[g]),
      seed = seed,
      quantization_error = mean(model$distances, na.rm = TRUE),
      occupied_nodes = length(unique(model$unit.classif))
    )
    for (k in K_VALUES) assignments[[paste0("k", k)]][, counter] <- partition_som(model, k)
  }
  summaries <- list(); sample <- list(); consensus <- list()
  for (k in K_VALUES) {
    key <- paste0("k", k)
    consensus[[key]] <- consensus_partition(assignments[[key]], k)
    labels <- consensus[[key]]$labels
    aris <- apply(assignments[[key]], 2, function(z) mclust::adjustedRandIndex(labels, z))
    sample_metrics <- sample_consensus_metrics(labels, consensus[[key]]$coassociation)
    summaries[[key]] <- data.frame(
      analysis = name,
      n = nrow(data),
      variables = paste(vars, collapse = "; "),
      k = k,
      models = ncol(assignments[[key]]),
      median_ari_to_consensus = median(aris),
      lower_ari = quantile(aris, 0.025),
      upper_ari = quantile(aris, 0.975),
      median_sample_support = median(sample_metrics$support),
      median_sample_entropy = median(sample_metrics$entropy)
    )
    sample[[key]] <- sample_metrics %>% mutate(analysis = name, k = k, .before = 1)
  }
  list(
    assignments = assignments,
    consensus = consensus,
    summary = bind_rows(summaries),
    sample_metrics = bind_rows(sample),
    model_metrics = bind_rows(metric_rows)
  )
}

module_results <- list(
  field4 = fit_som_sensitivity_ensemble(prepared$cohorts$field4, FIELD_VARS, character(), "field4", FALSE),
  nutrient4 = fit_som_sensitivity_ensemble(prepared$cohorts$nutrient4, NUTRIENT_VARS, NUTRIENT_VARS, "nutrient4", FALSE),
  nutrient_chla5 = fit_som_sensitivity_ensemble(prepared$cohorts$nutrient_chla5, NUTRIENT_CHLA_VARS, NUTRIENT_CHLA_VARS, "nutrient_chla5", FALSE),
  raw_core9 = fit_som_sensitivity_ensemble(prepared$cohorts$core9, VARS, character(), "raw_core9", FALSE),
  partial5 = fit_som_sensitivity_ensemble(prepared$cohorts$partial5, VARS, LOG_VARS, "partial5_native_na", TRUE)
)
message("Completed five SOM sensitivity ensembles")

data <- prepared$cohorts$core9
keep <- data$campaign != "2024-2025"
reduced <- data[keep, , drop = FALSE]
transform <- transform_fit(reduced)
x_reduced <- transform_apply(reduced, transform)
x_shared <- x_reduced

exclusion_rows <- list(); failures <- list(); counter <- 0L; failure_counter <- 0L
for (k in K_VALUES) {
  key <- paste0("k", k)
  methods <- list()
  set.seed(92000 + k)
  som <- kohonen::som(
    x_reduced,
    grid = kohonen::somgrid(8, 6, topo = "hexagonal"),
    rlen = SOM_RLEN, alpha = c(0.05, 0.01), keep.data = TRUE, mode = "online"
  )
  methods$SOM <- partition_som(som, k)
  methods[["K-means"]] <- fit_kmeans_predict(x_reduced, x_shared, k, 93000 + k)$labels
  methods[["Ward.D2"]] <- fit_ward_predict(x_reduced, x_shared, k)$labels
  gmm <- fit_gmm_predict(x_reduced, x_shared, k)
  if (is.na(gmm$failure)) methods$GMM <- gmm$labels else {
    failure_counter <- failure_counter + 1L
    failures[[failure_counter]] <- data.frame(analysis = "exclude_latest_campaign", method = "GMM", k = k, message = gmm$failure)
  }
  for (method in names(methods)) {
    reference <- resampling$consensus[[method]][[key]]$labels[keep]
    candidate <- methods[[method]]
    counter <- counter + 1L
    exclusion_rows[[counter]] <- data.frame(
      analysis = "exclude_latest_campaign",
      method = method,
      k = k,
      n = sum(keep),
      ari_to_primary_consensus = mclust::adjustedRandIndex(reference, candidate),
      ami_to_primary_consensus = adjusted_mutual_information(reference, candidate)
    )
  }
}

sensitivity <- list(
  modules = module_results,
  module_summary = bind_rows(lapply(module_results, `[[`, "summary")),
  exclusion = bind_rows(exclusion_rows),
  failures = if (length(failures)) bind_rows(failures) else data.frame(analysis = character(), method = character(), k = integer(), message = character())
)
save_rds_safe(sensitivity, file.path(MODEL_DIR, "sensitivity_ensembles.rds"))
write_csv_safe(sensitivity$module_summary, file.path(TABLE_DIR, "TableS_som_sensitivity_summary.csv"))
write_csv_safe(sensitivity$exclusion, file.path(TABLE_DIR, "TableS_latest_campaign_exclusion.csv"))
write_csv_safe(sensitivity$failures, file.path(TABLE_DIR, "TableS_sensitivity_failures.csv"))

cat("Sensitivity ensemble models:", sum(vapply(module_results, function(z) nrow(z$model_metrics), integer(1))), "\n")
cat("Exclusion comparison rows:", nrow(sensitivity$exclusion), "\n")
