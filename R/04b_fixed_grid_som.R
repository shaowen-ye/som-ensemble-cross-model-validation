args <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "00_config.R"))

# Like-for-like stability comparison.
#
# In 04_resampling_validation.R the SOM cycles through five grids across the 200
# resamples and takes its initialisation from the resampling seed, so its
# stability distribution absorbs grid and initialisation variation that the
# baselines never see. This script repeats the SOM arm on the identical
# resampling indices with the grid and the initialisation seed both held fixed,
# so that the only remaining source of variation is sample composition. The
# result is directly comparable with the K-means, Ward.D2 and GMM columns.

assert_packages(c("kohonen", "mclust", "cluster"))
prepared <- readRDS(file.path(RESULT_DIR, "prepared_data.rds"))
data <- prepared$cohorts$core9
n <- nrow(data)
n_resample <- as.integer(Sys.getenv("EI_N_RESAMPLE", unset = as.character(N_RESAMPLE)))
if (!is.finite(n_resample) || n_resample < 1L || n_resample > N_RESAMPLE) stop("Invalid EI_N_RESAMPLE")
seeds <- RESAMPLE_SEEDS[seq_len(n_resample)]

assignments <- setNames(
  lapply(K_VALUES, function(k) matrix(NA_integer_, n, n_resample)),
  paste0("k", K_VALUES)
)
failure_rows <- list()
failure_counter <- 0L

for (b in seq_len(n_resample)) {
  train_idx <- stratified_block_sample(data, RESAMPLE_FRACTION, seeds[b])
  train <- data[train_idx, , drop = FALSE]
  transform <- transform_fit(train)
  x_train <- transform_apply(train, transform)
  x_all <- transform_apply(data, transform)

  set.seed(FIXED_GRID_SEED)
  som_fit <- tryCatch(
    kohonen::som(
      x_train,
      grid = kohonen::somgrid(FIXED_GRID[["xdim"]], FIXED_GRID[["ydim"]], topo = "hexagonal"),
      rlen = SOM_RLEN,
      alpha = c(0.05, 0.01),
      keep.data = TRUE,
      mode = "online"
    ),
    error = function(e) e
  )

  for (k in K_VALUES) {
    key <- paste0("k", k)
    if (inherits(som_fit, "error")) {
      failure_counter <- failure_counter + 1L
      failure_rows[[failure_counter]] <- data.frame(
        resample = b, method = "SOM (fixed grid)", k = k, message = conditionMessage(som_fit)
      )
    } else {
      assignments[[key]][, b] <- partition_som(som_fit, k, x_all)
    }
  }
  if (b %% 20L == 0L || b == n_resample) message("Fixed-grid resampling: ", b, "/", n_resample)
}

consensus <- list()
stability_rows <- list()
cluster_rows <- list()
validity_rows <- list()

for (k in K_VALUES) {
  key <- paste0("k", k)
  cons <- consensus_partition(assignments[[key]], k)
  consensus[[key]] <- cons
  labels <- cons$labels

  replicate_stats <- lapply(seq_len(n_resample), function(b) {
    candidate <- assignments[[key]][, b]
    if (all(is.na(candidate))) return(NULL)
    jac <- clusterwise_jaccard(labels, candidate)
    data.frame(
      method = "SOM (fixed grid)",
      k = k,
      resample = b,
      ari = mclust::adjustedRandIndex(labels, candidate),
      ami = adjusted_mutual_information(labels, candidate),
      mean_cluster_jaccard = mean(jac$jaccard),
      minimum_cluster_jaccard = min(jac$jaccard)
    )
  })
  stability_rows[[key]] <- bind_rows(replicate_stats)

  cluster_rep <- lapply(seq_len(n_resample), function(b) {
    candidate <- assignments[[key]][, b]
    if (all(is.na(candidate))) return(NULL)
    clusterwise_jaccard(labels, candidate) %>% mutate(resample = b, .before = 1)
  }) %>% bind_rows()
  if (nrow(cluster_rep)) {
    cluster_rows[[key]] <- cluster_rep %>%
      group_by(cluster) %>%
      summarise(
        median_jaccard = median(jaccard),
        lower_95 = quantile(jaccard, 0.025),
        upper_95 = quantile(jaccard, 0.975),
        .groups = "drop"
      ) %>%
      mutate(method = "SOM (fixed grid)", k = k, stability = as.character(stability_class(median_jaccard)), .before = 1)
  }

  validity_rows[[key]] <- validity_metrics(prepared$x, labels) %>%
    mutate(method = "SOM (fixed grid)", k = k, .before = 1)
}

fixed_grid <- list(
  assignments = assignments,
  consensus = consensus,
  replicate_stability = bind_rows(stability_rows),
  cluster_stability = bind_rows(cluster_rows),
  internal_validity = bind_rows(validity_rows),
  failures = if (length(failure_rows)) bind_rows(failure_rows) else
    data.frame(resample = integer(), method = character(), k = integer(), message = character()),
  settings = list(
    grid = sprintf("%dx%d", FIXED_GRID[["xdim"]], FIXED_GRID[["ydim"]]),
    init_seed = FIXED_GRID_SEED,
    n_resample = n_resample,
    fraction = RESAMPLE_FRACTION,
    seeds = seeds
  )
)

suffix <- if (n_resample == N_RESAMPLE) "200" else paste0("test", n_resample)
save_rds_safe(fixed_grid, file.path(MODEL_DIR, paste0("fixed_grid_som_", suffix, ".rds")))
write_csv_safe(fixed_grid$replicate_stability, file.path(TABLE_DIR, paste0("TableS_fixed_grid_stability_", suffix, ".csv")))
write_csv_safe(fixed_grid$cluster_stability, file.path(TABLE_DIR, paste0("TableS_fixed_grid_clusterwise_jaccard_", suffix, ".csv")))
write_csv_safe(fixed_grid$internal_validity, file.path(TABLE_DIR, paste0("TableS_fixed_grid_internal_validity_", suffix, ".csv")))
write_csv_safe(fixed_grid$failures, file.path(TABLE_DIR, paste0("TableS_fixed_grid_failures_", suffix, ".csv")))

cat("Fixed-grid SOM:", fixed_grid$settings$grid, "over", n_resample, "identical grouped resamples\n")
cat("Recorded failures:", nrow(fixed_grid$failures), "\n")
