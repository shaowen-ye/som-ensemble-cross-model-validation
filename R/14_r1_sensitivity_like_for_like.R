args <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "00_config.R"))

# Like-for-like sensitivity comparisons (first-round revision).
#
# The SOM sensitivity ensembles of 06_sensitivity.R are grid-seed ensembles
# (five grids x ten seeds) fitted to complete cohorts, whereas the primary SOM
# consensus comes from 200 grouped resamples. This script supplies the matching
# grid-seed baseline for the nine-variable matrix, refits the variable modules on
# the same 390 events so that variables and observations are changed one at a
# time, and repeats the 2024-2025 exclusion with a grid-seed ensemble instead of
# the single map used in 06. Removing 16 randomly chosen events gives the
# reference against which that exclusion is judged.

assert_packages(c("kohonen", "mclust", "parallel"))
prepared <- readRDS(file.path(RESULT_DIR, "prepared_data.rds"))
ensemble <- readRDS(file.path(MODEL_DIR, "som_ensemble_250.rds"))
resampling <- readRDS(file.path(MODEL_DIR, "cross_model_resampling_200.rds"))
sensitivity <- readRDS(file.path(MODEL_DIR, "sensitivity_ensembles.rds"))
data <- prepared$cohorts$core9
methods <- c("SOM", "K-means", "Ward.D2", "GMM")
cores <- as.integer(Sys.getenv("EI_CORES", unset = as.character(max(1L, min(4L, parallel::detectCores() - 1L)))))

fit_grid_seed_ensemble <- function(x) {
  n_models <- nrow(GRIDS) * length(SENSITIVITY_SEEDS)
  assignments <- setNames(lapply(K_VALUES, function(k) matrix(NA_integer_, nrow(x), n_models)), paste0("k", K_VALUES))
  counter <- 0L
  for (g in seq_len(nrow(GRIDS))) for (seed in SENSITIVITY_SEEDS) {
    counter <- counter + 1L
    set.seed(seed)
    model <- kohonen::som(
      x,
      grid = kohonen::somgrid(GRIDS$xdim[g], GRIDS$ydim[g], topo = "hexagonal"),
      rlen = SOM_RLEN, alpha = c(0.05, 0.01), keep.data = TRUE, mode = "online"
    )
    for (k in K_VALUES) assignments[[paste0("k", k)]][, counter] <- partition_som(model, k)
  }
  assignments
}

ensemble_profile <- function(assignments, k) {
  labels <- consensus_partition(assignments, k)$labels
  aris <- apply(assignments, 2, function(z) mclust::adjustedRandIndex(labels, z))
  list(labels = labels, within = data.frame(
    within_median = median(aris),
    within_lower = unname(quantile(aris, 0.025)),
    within_upper = unname(quantile(aris, 0.975))
  ))
}

ari <- mclust::adjustedRandIndex

# Primary grid-seed baseline: the 50 maps of the 250-map ensemble that use the
# sensitivity seeds, in the same grid-major order as fit_grid_seed_ensemble().
baseline_columns <- which(ensemble$metrics$seed %in% SENSITIVITY_SEEDS)
if (length(baseline_columns) != nrow(GRIDS) * length(SENSITIVITY_SEEDS)) stop("Unexpected baseline map count")
baseline <- lapply(K_VALUES, function(k) ensemble_profile(ensemble$assignments[[paste0("k", k)]][, baseline_columns, drop = FALSE], k))
names(baseline) <- paste0("k", K_VALUES)
baseline_table <- bind_rows(lapply(K_VALUES, function(k) {
  key <- paste0("k", k)
  baseline[[key]]$within %>% mutate(
    analysis = "Primary nine-variable matrix", cohort = "Primary 390 events", n = nrow(data), k = k,
    agreement_primary_grid_seed = 1,
    agreement_primary_resampling = ari(baseline[[key]]$labels, resampling$consensus$SOM[[key]]$labels),
    .before = 1
  )
}))

# Variable modules on the common 390 events, and the original module cohorts ----
module_specs <- list(
  field4 = list(label = "Field chemistry", vars = FIELD_VARS, log_vars = character()),
  nutrient4 = list(label = "Nutrients", vars = NUTRIENT_VARS, log_vars = NUTRIENT_VARS),
  nutrient_chla5 = list(label = "Chl-a + nutrients", vars = NUTRIENT_CHLA_VARS, log_vars = NUTRIENT_CHLA_VARS)
)
cohort_data <- list(field4 = prepared$cohorts$field4, nutrient4 = prepared$cohorts$nutrient4,
                    nutrient_chla5 = prepared$cohorts$nutrient_chla5, raw_core9 = prepared$cohorts$core9,
                    partial5 = prepared$cohorts$partial5)
module_rows <- list(); common_consensus <- list(); counter <- 0L
for (name in names(module_specs)) {
  spec <- module_specs[[name]]
  x <- transform_apply(data, transform_fit(data, vars = spec$vars, log_vars = spec$log_vars))
  assignments <- fit_grid_seed_ensemble(x)
  index <- match(data$record_id, cohort_data[[name]]$record_id)
  for (k in K_VALUES) {
    key <- paste0("k", k)
    profile <- ensemble_profile(assignments[[key]], k)
    common_consensus[[paste(name, key)]] <- profile$labels
    full_labels <- sensitivity$modules[[name]]$consensus[[key]]$labels[index]
    full_summary <- sensitivity$module_summary %>% filter(analysis == name, k == !!k)
    counter <- counter + 1L
    module_rows[[counter]] <- profile$within %>% mutate(
      analysis = spec$label, cohort = "Primary 390 events", n = nrow(data), k = k,
      agreement_primary_grid_seed = ari(profile$labels, baseline[[key]]$labels),
      agreement_primary_resampling = ari(profile$labels, resampling$consensus$SOM[[key]]$labels),
      agreement_with_full_module_cohort = ari(profile$labels, full_labels),
      .before = 1
    )
    counter <- counter + 1L
    module_rows[[counter]] <- data.frame(
      analysis = spec$label, cohort = "Full module cohort", n = nrow(cohort_data[[name]]), k = k,
      within_median = full_summary$median_ari_to_consensus, within_lower = full_summary$lower_ari, within_upper = full_summary$upper_ari,
      agreement_primary_grid_seed = ari(full_labels, baseline[[key]]$labels),
      agreement_primary_resampling = ari(full_labels, resampling$consensus$SOM[[key]]$labels),
      agreement_with_full_module_cohort = 1
    )
  }
  message("Module on common events completed: ", name)
}
for (name in c("raw_core9", "partial5")) {
  label <- c(raw_core9 = "Untransformed nine variables", partial5 = "Native missing distance")[[name]]
  summary_name <- c(raw_core9 = "raw_core9", partial5 = "partial5_native_na")[[name]]
  index <- match(data$record_id, cohort_data[[name]]$record_id)
  for (k in K_VALUES) {
    key <- paste0("k", k)
    full_labels <- sensitivity$modules[[name]]$consensus[[key]]$labels[index]
    full_summary <- sensitivity$module_summary %>% filter(analysis == summary_name, k == !!k)
    counter <- counter + 1L
    module_rows[[counter]] <- data.frame(
      analysis = label, cohort = if (name == "partial5") "Partial cohort (454 events)" else "Primary 390 events",
      n = nrow(cohort_data[[name]]), k = k,
      within_median = full_summary$median_ari_to_consensus, within_lower = full_summary$lower_ari, within_upper = full_summary$upper_ari,
      agreement_primary_grid_seed = ari(full_labels, baseline[[key]]$labels),
      agreement_primary_resampling = ari(full_labels, resampling$consensus$SOM[[key]]$labels),
      agreement_with_full_module_cohort = NA_real_
    )
  }
}
module_table <- bind_rows(baseline_table %>% mutate(agreement_with_full_module_cohort = NA_real_), bind_rows(module_rows))

# The native missing-distance fit reduces to the primary SOM on complete events:
# the initialisation draw and the online updates use the same random stream.
# Checked for all 50 maps of the sensitivity design (five grids x ten seeds).
x_primary <- prepared$x
identity_rows <- list()
for (g in seq_len(nrow(GRIDS))) for (seed in SENSITIVITY_SEEDS) {
  nodes <- GRIDS$xdim[g] * GRIDS$ydim[g]
  set.seed(seed)
  init <- x_primary[sample(seq_len(nrow(x_primary)), nodes, replace = FALSE), , drop = FALSE]
  model <- kohonen::som(x_primary, grid = kohonen::somgrid(GRIDS$xdim[g], GRIDS$ydim[g], topo = "hexagonal"),
                        rlen = SOM_RLEN, alpha = c(0.05, 0.01), keep.data = TRUE, mode = "online",
                        maxNA.fraction = 4 / 9, init = init)
  column <- which(ensemble$metrics$grid == sprintf("%dx%d", GRIDS$xdim[g], GRIDS$ydim[g]) & ensemble$metrics$seed == seed)
  identity_rows[[length(identity_rows) + 1L]] <- data.frame(
    grid = sprintf("%dx%d", GRIDS$xdim[g], GRIDS$ydim[g]), seed = seed,
    identical_partitions_k2_8 = all(vapply(K_VALUES, function(k) identical(
      as.integer(partition_som(model, k)), as.integer(ensemble$assignments[[paste0("k", k)]][, column])), logical(1))),
    max_codebook_difference = max(abs(model$codes[[1]] - ensemble$models[[column]]$codes[[1]])))
}
missing_identity <- bind_rows(identity_rows)
if (!all(missing_identity$identical_partitions_k2_8) || max(missing_identity$max_codebook_difference) != 0) {
  stop("The missing-distance SOM does not reproduce the primary maps on complete events")
}

# Exclusion of the 2024-2025 survey with a grid-seed ensemble -------------------
keep <- data$campaign != "2024-2025"
reduced <- data[keep, , drop = FALSE]
x_reduced <- transform_apply(reduced, transform_fit(reduced))
reduced_assignments <- fit_grid_seed_ensemble(x_reduced)
exclusion_rows <- list(); counter <- 0L
for (k in K_VALUES) {
  key <- paste0("k", k)
  profile <- ensemble_profile(reduced_assignments[[key]], k)
  reference <- resampling$consensus$SOM[[key]]$labels[keep]
  singles <- apply(reduced_assignments[[key]], 2, function(z) ari(reference, z))
  original <- sensitivity$exclusion %>% filter(method == "SOM", k == !!k)
  counter <- counter + 1L
  exclusion_rows[[counter]] <- data.frame(
    method = "SOM", k = k, n = sum(keep),
    submitted_single_map_ari = original$ari_to_primary_consensus,
    ensemble_consensus_ari = ari(reference, profile$labels),
    ensemble_vs_primary_grid_seed = ari(baseline[[key]]$labels[keep], profile$labels),
    single_map_median = median(singles), single_map_lower = unname(quantile(singles, 0.025)),
    single_map_upper = unname(quantile(singles, 0.975))
  )
  for (method in setdiff(methods, "SOM")) {
    original <- sensitivity$exclusion %>% filter(method == !!method, k == !!k)
    counter <- counter + 1L
    exclusion_rows[[counter]] <- data.frame(
      method = method, k = k, n = sum(keep),
      submitted_single_map_ari = NA_real_, ensemble_consensus_ari = original$ari_to_primary_consensus,
      ensemble_vs_primary_grid_seed = NA_real_, single_map_median = NA_real_, single_map_lower = NA_real_, single_map_upper = NA_real_
    )
  }
}
exclusion <- bind_rows(exclusion_rows)

# Reference: remove randomly chosen blocks totalling about 16 events -------------
target <- sum(!keep)
run_omission <- function(r) {
  set.seed(EXCLUSION_OMISSION_SEEDS[r])
  # Blocks are listed in order of first appearance, so the draw does not depend
  # on how block identifiers are spelled or sorted in the current locale.
  blocks <- data.frame(block_id = unique(data$block_id), stringsAsFactors = FALSE)
  blocks$events <- as.integer(table(data$block_id)[blocks$block_id])
  blocks <- blocks[sample.int(nrow(blocks)), , drop = FALSE]
  cut <- which.min(abs(cumsum(blocks$events) - target))
  retained <- !(data$block_id %in% blocks$block_id[seq_len(cut)])
  sub <- data[retained, , drop = FALSE]
  x_sub <- transform_apply(sub, transform_fit(sub))
  som_assignments <- fit_grid_seed_ensemble(x_sub)
  rows <- lapply(K_VALUES, function(k) {
    key <- paste0("k", k)
    gmm <- fit_gmm_predict(x_sub, x_sub, k)
    data.frame(
      replicate = r, omitted_n = sum(!retained), k = k,
      SOM = ari(resampling$consensus$SOM[[key]]$labels[retained], consensus_partition(som_assignments[[key]], k)$labels),
      `K-means` = ari(resampling$consensus[["K-means"]][[key]]$labels[retained], fit_kmeans_predict(x_sub, x_sub, k, 93000 + k)$labels),
      Ward.D2 = ari(resampling$consensus[["Ward.D2"]][[key]]$labels[retained], fit_ward_predict(x_sub, x_sub, k)$labels),
      GMM = if (is.na(gmm$failure)) ari(resampling$consensus$GMM[[key]]$labels[retained], gmm$labels) else NA_real_,
      check.names = FALSE
    )
  })
  bind_rows(rows)
}
omission <- parallel::mclapply(seq_len(MATCHED_OMISSION_REPLICATES), run_omission, mc.cores = cores, mc.preschedule = FALSE)
if (any(vapply(omission, inherits, logical(1), what = "try-error"))) stop("Random-omission replicates failed")
omission <- bind_rows(omission) %>%
  tidyr::pivot_longer(all_of(methods), names_to = "method", values_to = "ari")
omission_summary <- omission %>%
  group_by(method, k) %>%
  summarise(random_omission_median = median(ari, na.rm = TRUE),
            random_omission_lower = unname(quantile(ari, 0.025, na.rm = TRUE)),
            random_omission_upper = unname(quantile(ari, 0.975, na.rm = TRUE)),
            omitted_n_median = median(omitted_n), .groups = "drop")
exclusion <- exclusion %>% left_join(omission_summary, by = c("method", "k"))

like_for_like <- list(
  baseline = baseline_table,
  modules = module_table,
  common_consensus = common_consensus,
  missing_identity = missing_identity,
  exclusion = exclusion,
  exclusion_omission = omission
)
save_rds_safe(like_for_like, file.path(MODEL_DIR, "r1_sensitivity_like_for_like.rds"))
write_csv_safe(module_table, file.path(TABLE_DIR, "TableS_r1_modules_like_for_like.csv"))
write_csv_safe(missing_identity, file.path(TABLE_DIR, "TableS_r1_missing_distance_identity_check.csv"))
write_csv_safe(exclusion, file.path(TABLE_DIR, "TableS_r1_survey_exclusion.csv"))

cat("Like-for-like sensitivity completed\n")
print(module_table %>% filter(k == 3L) %>% select(analysis, cohort, n, within_median, agreement_primary_grid_seed, agreement_primary_resampling, agreement_with_full_module_cohort))
print(missing_identity)
print(exclusion %>% filter(k %in% 2:4))
