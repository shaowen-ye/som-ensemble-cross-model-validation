args <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "00_config.R"))

# Heavy-tail sensitivity.
#
# The primary matrix retains genuinely extreme observations, most of them low
# nitrate values recorded during summer algal blooms. Because SOM prototypes and
# Ward.D2 tree cuts can isolate such points into very small clusters, this script
# repeats the full four-model grouped resampling after removing every event whose
# standardised value on any variable exceeds OUTLIER_Z_THRESHOLD. It tests
# directly whether the fragile-cluster result is driven by heavy tails rather
# than by the structure of the data.

assert_packages(c("kohonen", "mclust", "cluster"))
prepared <- readRDS(file.path(RESULT_DIR, "prepared_data.rds"))
full_data <- prepared$cohorts$core9
keep <- apply(abs(prepared$x), 1, max) <= OUTLIER_Z_THRESHOLD
data <- full_data[keep, , drop = FALSE]
n <- nrow(data)
dropped <- full_data[!keep, c("sampling_event_id", "campaign"), drop = FALSE]
dropped$max_abs_z <- round(apply(abs(prepared$x[!keep, , drop = FALSE]), 1, max), 3)
dropped$driving_variable <- colnames(prepared$x)[apply(abs(prepared$x[!keep, , drop = FALSE]), 1, which.max)]
write_csv_safe(dropped, file.path(TABLE_DIR, "TableS_outlier_dropped_events.csv"))

n_resample <- as.integer(Sys.getenv("EI_N_RESAMPLE", unset = as.character(N_RESAMPLE)))
if (!is.finite(n_resample) || n_resample < 1L || n_resample > N_RESAMPLE) stop("Invalid EI_N_RESAMPLE")
seeds <- RESAMPLE_SEEDS[seq_len(n_resample)]
methods <- c("SOM", "K-means", "Ward.D2", "GMM")

assignments <- setNames(lapply(methods, function(method) {
  setNames(lapply(K_VALUES, function(k) matrix(NA_integer_, n, n_resample)), paste0("k", K_VALUES))
}), methods)
failure_rows <- list()
failure_counter <- 0L

for (b in seq_len(n_resample)) {
  train_idx <- stratified_block_sample(data, RESAMPLE_FRACTION, seeds[b])
  train <- data[train_idx, , drop = FALSE]
  transform <- transform_fit(train)
  x_train <- transform_apply(train, transform)
  x_all <- transform_apply(data, transform)

  g <- 1L + ((b - 1L) %% nrow(GRIDS))
  set.seed(seeds[b])
  som_fit <- tryCatch(
    kohonen::som(
      x_train,
      grid = kohonen::somgrid(GRIDS$xdim[g], GRIDS$ydim[g], topo = "hexagonal"),
      rlen = SOM_RLEN, alpha = c(0.05, 0.01), keep.data = TRUE, mode = "online"
    ),
    error = function(e) e
  )

  for (k in K_VALUES) {
    key <- paste0("k", k)
    if (inherits(som_fit, "error")) {
      failure_counter <- failure_counter + 1L
      failure_rows[[failure_counter]] <- data.frame(resample = b, method = "SOM", k = k, message = conditionMessage(som_fit))
    } else {
      assignments$SOM[[key]][, b] <- partition_som(som_fit, k, x_all)
    }

    km <- tryCatch(fit_kmeans_predict(x_train, x_all, k, seeds[b] + k),
                   error = function(e) list(labels = rep(NA_integer_, n), failure = conditionMessage(e)))
    assignments[["K-means"]][[key]][, b] <- km$labels
    if (!is.na(km$failure)) {
      failure_counter <- failure_counter + 1L
      failure_rows[[failure_counter]] <- data.frame(resample = b, method = "K-means", k = k, message = km$failure)
    }

    ward <- tryCatch(fit_ward_predict(x_train, x_all, k),
                     error = function(e) list(labels = rep(NA_integer_, n), failure = conditionMessage(e)))
    assignments[["Ward.D2"]][[key]][, b] <- ward$labels
    if (!is.na(ward$failure)) {
      failure_counter <- failure_counter + 1L
      failure_rows[[failure_counter]] <- data.frame(resample = b, method = "Ward.D2", k = k, message = ward$failure)
    }

    gmm <- fit_gmm_predict(x_train, x_all, k)
    assignments$GMM[[key]][, b] <- gmm$labels
    if (!is.na(gmm$failure)) {
      failure_counter <- failure_counter + 1L
      failure_rows[[failure_counter]] <- data.frame(resample = b, method = "GMM", k = k, message = gmm$failure)
    }
  }
  if (b %% 20L == 0L || b == n_resample) message("Outlier-trimmed resampling: ", b, "/", n_resample)
}

stability_rows <- list()
cluster_rows <- list()

for (method in methods) {
  for (k in K_VALUES) {
    key <- paste0("k", k)
    cons <- consensus_partition(assignments[[method]][[key]], k)
    labels <- cons$labels

    replicate_stats <- lapply(seq_len(n_resample), function(b) {
      candidate <- assignments[[method]][[key]][, b]
      if (all(is.na(candidate))) return(NULL)
      jac <- clusterwise_jaccard(labels, candidate)
      data.frame(
        method = method, k = k, resample = b,
        ari = mclust::adjustedRandIndex(labels, candidate),
        ami = adjusted_mutual_information(labels, candidate),
        mean_cluster_jaccard = mean(jac$jaccard),
        minimum_cluster_jaccard = min(jac$jaccard)
      )
    })
    stability_rows[[paste(method, k)]] <- bind_rows(replicate_stats)

    cluster_rep <- lapply(seq_len(n_resample), function(b) {
      candidate <- assignments[[method]][[key]][, b]
      if (all(is.na(candidate))) return(NULL)
      clusterwise_jaccard(labels, candidate) %>% mutate(resample = b, .before = 1)
    }) %>% bind_rows()
    if (nrow(cluster_rep)) {
      cluster_rows[[paste(method, k)]] <- cluster_rep %>%
        group_by(cluster) %>%
        summarise(
          median_jaccard = median(jaccard),
          lower_95 = quantile(jaccard, 0.025),
          upper_95 = quantile(jaccard, 0.975),
          .groups = "drop"
        ) %>%
        mutate(method = method, k = k, stability = as.character(stability_class(median_jaccard)), .before = 1)
    }
  }
}

outlier_sensitivity <- list(
  replicate_stability = bind_rows(stability_rows),
  cluster_stability = bind_rows(cluster_rows),
  dropped_events = dropped,
  failures = if (length(failure_rows)) bind_rows(failure_rows) else
    data.frame(resample = integer(), method = character(), k = integer(), message = character()),
  settings = list(z_threshold = OUTLIER_Z_THRESHOLD, n_retained = n, n_dropped = nrow(dropped), n_resample = n_resample)
)

suffix <- if (n_resample == N_RESAMPLE) "200" else paste0("test", n_resample)
save_rds_safe(outlier_sensitivity, file.path(MODEL_DIR, paste0("outlier_sensitivity_", suffix, ".rds")))
write_csv_safe(outlier_sensitivity$replicate_stability, file.path(TABLE_DIR, paste0("TableS_outlier_sensitivity_stability_", suffix, ".csv")))
write_csv_safe(outlier_sensitivity$cluster_stability, file.path(TABLE_DIR, paste0("TableS_outlier_sensitivity_clusterwise_jaccard_", suffix, ".csv")))
write_csv_safe(outlier_sensitivity$failures, file.path(TABLE_DIR, paste0("TableS_outlier_sensitivity_failures_", suffix, ".csv")))

cat("Outlier-trimmed cohort:", n, "of", nrow(full_data), "events retained at |z| <=", OUTLIER_Z_THRESHOLD, "\n")
cat("Recorded failures:", nrow(outlier_sensitivity$failures), "\n")
