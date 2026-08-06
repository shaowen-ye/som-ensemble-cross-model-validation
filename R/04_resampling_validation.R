# Compare SOM, K-means, Ward.D2 and GMM under identical grouped resamples.
args <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "00_config.R"))

assert_packages(c("kohonen", "mclust", "cluster"))
prepared <- readRDS(file.path(RESULT_DIR, "prepared_data.rds"))
data <- prepared$cohorts$core9
n <- nrow(data)
n_resample <- as.integer(Sys.getenv("EI_N_RESAMPLE", unset = as.character(N_RESAMPLE)))
if (!is.finite(n_resample) || n_resample < 1L || n_resample > N_RESAMPLE) stop("Invalid EI_N_RESAMPLE")
seeds <- RESAMPLE_SEEDS[seq_len(n_resample)]
methods <- c("SOM", "K-means", "Ward.D2", "GMM")

# Store full-sample predictions from every resampling replicate and method.
assignments <- setNames(lapply(methods, function(method) {
  setNames(lapply(K_VALUES, function(k) matrix(NA_integer_, n, n_resample)), paste0("k", K_VALUES))
}), methods)

resample_manifest <- vector("list", n_resample)
failure_rows <- list()
gmm_rows <- list()
failure_counter <- 0L
gmm_counter <- 0L

# Draw 80% of campaign-specific station blocks and refit all four model families.
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
      rlen = SOM_RLEN,
      alpha = c(0.05, 0.01),
      keep.data = TRUE,
      mode = "online"
    ),
    error = function(e) e
  )

  resample_manifest[[b]] <- data.frame(
    resample = b,
    seed = seeds[b],
    n_train = length(train_idx),
    blocks_train = n_distinct(train$block_id),
    som_grid = sprintf("%dx%d", GRIDS$xdim[g], GRIDS$ydim[g])
  )

  for (k in K_VALUES) {
    key <- paste0("k", k)
    if (inherits(som_fit, "error")) {
      failure_counter <- failure_counter + 1L
      failure_rows[[failure_counter]] <- data.frame(resample = b, method = "SOM", k = k, message = conditionMessage(som_fit))
    } else {
      assignments$SOM[[key]][, b] <- partition_som(som_fit, k, x_all)
    }

    km <- tryCatch(fit_kmeans_predict(x_train, x_all, k, seeds[b] + k), error = function(e) list(labels = rep(NA_integer_, n), failure = conditionMessage(e)))
    assignments[["K-means"]][[key]][, b] <- km$labels
    if (!is.na(km$failure)) {
      failure_counter <- failure_counter + 1L
      failure_rows[[failure_counter]] <- data.frame(resample = b, method = "K-means", k = k, message = km$failure)
    }

    ward <- tryCatch(fit_ward_predict(x_train, x_all, k), error = function(e) list(labels = rep(NA_integer_, n), failure = conditionMessage(e)))
    assignments[["Ward.D2"]][[key]][, b] <- ward$labels
    if (!is.na(ward$failure)) {
      failure_counter <- failure_counter + 1L
      failure_rows[[failure_counter]] <- data.frame(resample = b, method = "Ward.D2", k = k, message = ward$failure)
    }

    gmm <- fit_gmm_predict(x_train, x_all, k)
    assignments$GMM[[key]][, b] <- gmm$labels
    gmm_counter <- gmm_counter + 1L
    gmm_rows[[gmm_counter]] <- data.frame(resample = b, k = k, model_name = gmm$model_name, converged = is.na(gmm$failure))
    if (!is.na(gmm$failure)) {
      failure_counter <- failure_counter + 1L
      failure_rows[[failure_counter]] <- data.frame(resample = b, method = "GMM", k = k, message = gmm$failure)
    }
  }
  if (b %% 5L == 0L || b == n_resample) message("Grouped resampling: ", b, "/", n_resample)
}

# Derive method-specific consensus partitions, stability and sample uncertainty.
consensus <- setNames(vector("list", length(methods)), methods)
sample_rows <- list()
stability_rows <- list()
cluster_rows <- list()
validity_rows <- list()
counter <- 0L

for (method in methods) {
  consensus[[method]] <- list()
  for (k in K_VALUES) {
    key <- paste0("k", k)
    cons <- consensus_partition(assignments[[method]][[key]], k)
    consensus[[method]][[key]] <- cons
    labels <- cons$labels
    sm <- sample_consensus_metrics(labels, cons$coassociation) %>%
      mutate(method = method, k = k, consensus_cluster = labels, .before = 1)
    sample_rows[[paste(method, k)]] <- sm

    replicate_stats <- lapply(seq_len(n_resample), function(b) {
      candidate <- assignments[[method]][[key]][, b]
      if (all(is.na(candidate))) return(NULL)
      jac <- clusterwise_jaccard(labels, candidate)
      data.frame(
        method = method,
        k = k,
        resample = b,
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

    counter <- counter + 1L
    validity_rows[[counter]] <- validity_metrics(prepared$x, labels) %>% mutate(method = method, k = k, .before = 1)
  }
}

# Quantify agreement across model families without constructing a ranking score.
cross_model_rows <- list()
counter <- 0L
for (k in K_VALUES) {
  pairs <- combn(methods, 2, simplify = FALSE)
  for (pair in pairs) {
    counter <- counter + 1L
    a <- consensus[[pair[1]]][[paste0("k", k)]]$labels
    b <- consensus[[pair[2]]][[paste0("k", k)]]$labels
    cross_model_rows[[counter]] <- data.frame(
      k = k,
      method_1 = pair[1],
      method_2 = pair[2],
      ari = mclust::adjustedRandIndex(a, b),
      ami = adjusted_mutual_information(a, b)
    )
  }
}

resampling <- list(
  assignments = assignments,
  consensus = consensus,
  sample_metrics = bind_rows(sample_rows),
  replicate_stability = bind_rows(stability_rows),
  cluster_stability = bind_rows(cluster_rows),
  cross_model_agreement = bind_rows(cross_model_rows),
  internal_validity = bind_rows(validity_rows),
  resample_manifest = bind_rows(resample_manifest),
  failures = if (length(failure_rows)) bind_rows(failure_rows) else data.frame(resample = integer(), method = character(), k = integer(), message = character()),
  gmm_models = bind_rows(gmm_rows),
  settings = list(n_resample = n_resample, fraction = RESAMPLE_FRACTION, seeds = seeds)
)

# Preserve failures and GMM covariance selections instead of silently dropping them.
suffix <- if (n_resample == N_RESAMPLE) "200" else paste0("test", n_resample)
save_rds_safe(resampling, file.path(MODEL_DIR, paste0("cross_model_resampling_", suffix, ".rds")))
write_csv_safe(resampling$resample_manifest, file.path(TABLE_DIR, paste0("TableS_resample_manifest_", suffix, ".csv")))
write_csv_safe(resampling$replicate_stability, file.path(TABLE_DIR, paste0("TableS_resampling_stability_", suffix, ".csv")))
write_csv_safe(resampling$cluster_stability, file.path(TABLE_DIR, paste0("TableS_clusterwise_jaccard_", suffix, ".csv")))
write_csv_safe(resampling$cross_model_agreement, file.path(TABLE_DIR, paste0("TableS_cross_model_agreement_", suffix, ".csv")))
write_csv_safe(resampling$internal_validity, file.path(TABLE_DIR, paste0("TableS_internal_validity_", suffix, ".csv")))
write_csv_safe(resampling$sample_metrics, file.path(TABLE_DIR, paste0("TableS_sample_consensus_metrics_", suffix, ".csv")))
write_csv_safe(resampling$failures, file.path(TABLE_DIR, paste0("TableS_model_failures_", suffix, ".csv")))
write_csv_safe(resampling$gmm_models, file.path(TABLE_DIR, paste0("TableS_gmm_models_", suffix, ".csv")))

cat("Completed", n_resample, "grouped resamples for four model families\n")
cat("Recorded failures:", nrow(resampling$failures), "\n")
