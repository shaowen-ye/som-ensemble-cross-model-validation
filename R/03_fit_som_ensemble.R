args <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "00_config.R"))

assert_packages(c("kohonen", "mclust", "cluster"))
prepared_file <- file.path(RESULT_DIR, "prepared_data.rds")
if (!file.exists(prepared_file)) stop("Run 02_prepare_data.R first")
prepared <- readRDS(prepared_file)
x <- prepared$x
n_models <- nrow(GRIDS) * length(SOM_SEEDS)

metrics <- vector("list", n_models)
models <- vector("list", n_models)
assignments <- setNames(lapply(K_VALUES, function(k) matrix(NA_integer_, nrow(x), n_models)), paste0("k", K_VALUES))
model_names <- character(n_models)

counter <- 0L
for (g in seq_len(nrow(GRIDS))) {
  for (seed in SOM_SEEDS) {
    counter <- counter + 1L
    set.seed(seed)
    model <- kohonen::som(
      x,
      grid = kohonen::somgrid(GRIDS$xdim[g], GRIDS$ydim[g], topo = "hexagonal"),
      rlen = SOM_RLEN,
      alpha = c(0.05, 0.01),
      keep.data = TRUE,
      mode = "online"
    )
    label <- sprintf("%02dx%02d_seed%04d", GRIDS$xdim[g], GRIDS$ydim[g], seed)
    model_names[counter] <- label
    models[[counter]] <- model
    metrics[[counter]] <- data.frame(
      model = label,
      grid = sprintf("%dx%d", GRIDS$xdim[g], GRIDS$ydim[g]),
      nodes = GRIDS$xdim[g] * GRIDS$ydim[g],
      seed = seed,
      quantization_error = mean(model$distances),
      topographic_error = topographic_error(model, x),
      neighbour_preservation_10 = neighbour_preservation(model, x, 10L),
      occupied_nodes = length(unique(model$unit.classif)),
      empty_nodes = GRIDS$xdim[g] * GRIDS$ydim[g] - length(unique(model$unit.classif))
    )
    for (k in K_VALUES) assignments[[paste0("k", k)]][, counter] <- partition_som(model, k)
    if (counter %% 10L == 0L) message("SOM ensemble: ", counter, "/", n_models)
  }
}

metrics <- bind_rows(metrics)
for (key in names(assignments)) colnames(assignments[[key]]) <- model_names

consensus <- list()
sample_rows <- list()
model_agreement <- list()
pairwise_rows <- list()
for (k in K_VALUES) {
  key <- paste0("k", k)
  consensus[[key]] <- consensus_partition(assignments[[key]], k)
  labels <- consensus[[key]]$labels
  sample_rows[[key]] <- sample_consensus_metrics(labels, consensus[[key]]$coassociation) %>%
    mutate(k = k, consensus_cluster = labels, .before = 1)
  agreement <- lapply(seq_len(ncol(assignments[[key]])), function(j) {
    candidate <- assignments[[key]][, j]
    jac <- clusterwise_jaccard(labels, candidate)
    data.frame(
      k = k,
      model = colnames(assignments[[key]])[j],
      ari_to_consensus = mclust::adjustedRandIndex(labels, candidate),
      ami_to_consensus = adjusted_mutual_information(labels, candidate),
      mean_cluster_jaccard = mean(jac$jaccard),
      minimum_cluster_jaccard = min(jac$jaccard)
    )
  })
  model_agreement[[key]] <- bind_rows(agreement)

  cmb <- combn(seq_len(ncol(assignments[[key]])), 2)
  pairwise_rows[[key]] <- data.frame(
    k = k,
    model_1 = colnames(assignments[[key]])[cmb[1, ]],
    model_2 = colnames(assignments[[key]])[cmb[2, ]],
    ari = apply(cmb, 2, function(ij) mclust::adjustedRandIndex(assignments[[key]][, ij[1]], assignments[[key]][, ij[2]]))
  )
}

model_agreement <- bind_rows(model_agreement)
metrics <- metrics %>%
  mutate(
    rank_qe = rank(quantization_error, ties.method = "average"),
    rank_te = rank(topographic_error, ties.method = "average"),
    rank_np = rank(-neighbour_preservation_10, ties.method = "average"),
    topology_rank_sum = rank_qe + rank_te + rank_np
  )
representative_index <- which.min(metrics$topology_rank_sum)

ensemble <- list(
  metrics = metrics,
  assignments = assignments,
  consensus = consensus,
  sample_metrics = bind_rows(sample_rows),
  model_agreement = model_agreement,
  pairwise_ari = bind_rows(pairwise_rows),
  models = models,
  representative_index = representative_index,
  representative_name = model_names[representative_index],
  settings = list(grids = GRIDS, seeds = SOM_SEEDS, rlen = SOM_RLEN, k = K_VALUES)
)

save_rds_safe(ensemble, file.path(MODEL_DIR, "som_ensemble_250.rds"))
write_csv_safe(metrics, file.path(TABLE_DIR, "TableS_som_ensemble_metrics.csv"))
write_csv_safe(model_agreement, file.path(TABLE_DIR, "TableS_som_model_consensus_agreement.csv"))
write_csv_safe(ensemble$pairwise_ari, file.path(TABLE_DIR, "TableS_som_pairwise_ari.csv"))
write_csv_safe(ensemble$sample_metrics, file.path(TABLE_DIR, "TableS_som_sample_support.csv"))

cat("Completed", n_models, "SOMs\n")
cat("Representative topology model:", ensemble$representative_name, "\n")
