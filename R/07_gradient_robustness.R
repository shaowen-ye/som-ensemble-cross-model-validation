args <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "00_config.R"))

prepared <- readRDS(file.path(RESULT_DIR, "prepared_data.rds"))
ensemble <- readRDS(file.path(MODEL_DIR, "som_ensemble_250.rds"))
data <- prepared$cohorts$core9
x_full <- prepared$x

full_pca <- prcomp(x_full, center = FALSE, scale. = FALSE)
reference <- scale(full_pca$x[, 1:2], center = TRUE, scale = FALSE)
reference_dist <- as.numeric(dist(reference))

align_scores <- function(target, reference) {
  target <- scale(target, center = TRUE, scale = FALSE)
  cross <- crossprod(target, reference)
  decomposition <- svd(cross)
  rotation <- decomposition$u %*% t(decomposition$v)
  rotated <- target %*% rotation
  scale_factor <- sum(reference * rotated) / sum(rotated^2)
  aligned <- rotated * scale_factor
  list(
    scores = aligned,
    rotation = rotation,
    scale_factor = scale_factor,
    correlation = sum(reference * aligned) / sqrt(sum(reference^2) * sum(aligned^2))
  )
}

pca_rows <- vector("list", N_RESAMPLE)
loading_rows <- vector("list", N_RESAMPLE)
for (b in seq_len(N_RESAMPLE)) {
  train_idx <- stratified_block_sample(data, RESAMPLE_FRACTION, RESAMPLE_SEEDS[b])
  train <- data[train_idx, , drop = FALSE]
  fitted_transform <- transform_fit(train)
  x_train <- transform_apply(train, fitted_transform)
  x_all <- transform_apply(data, fitted_transform)
  fit <- prcomp(x_train, center = FALSE, scale. = FALSE)
  raw_scores <- x_all %*% fit$rotation[, 1:2, drop = FALSE]
  aligned <- align_scores(raw_scores, reference)

  rotated_loadings <- fit$rotation[, 1:2, drop = FALSE] %*% aligned$rotation * aligned$scale_factor
  loading_similarity <- vapply(seq_len(nrow(rotated_loadings)), function(j) {
    cosine <- sum(rotated_loadings[j, ] * full_pca$rotation[j, 1:2]) /
      sqrt(sum(rotated_loadings[j, ]^2) * sum(full_pca$rotation[j, 1:2]^2))
    cosine
  }, numeric(1))

  pca_rows[[b]] <- data.frame(
    resample = b,
    seed = RESAMPLE_SEEDS[b],
    n_train = length(train_idx),
    procrustes_correlation = aligned$correlation,
    pairwise_distance_spearman = cor(reference_dist, as.numeric(dist(aligned$scores)), method = "spearman"),
    variance_pc1_pc2 = sum(fit$sdev[1:2]^2) / sum(fit$sdev^2)
  )
  loading_rows[[b]] <- data.frame(
    resample = b,
    variable = rownames(rotated_loadings),
    loading_vector_cosine = loading_similarity
  )
}

input_dist <- as.numeric(dist(x_full))
input_relation <- cor(x_full)
relation_index <- upper.tri(input_relation)
som_rows <- vector("list", length(ensemble$models))
for (i in seq_along(ensemble$models)) {
  model <- ensemble$models[[i]]
  grid_dist <- kohonen::unit.distances(model$grid)
  mapped_dist <- grid_dist[model$unit.classif, model$unit.classif]
  code_relation <- cor(model$codes[[1]])
  som_rows[[i]] <- data.frame(
    model = ensemble$metrics$model[i],
    grid = ensemble$metrics$grid[i],
    seed = ensemble$metrics$seed[i],
    input_to_map_distance_spearman = cor(input_dist, as.numeric(as.dist(mapped_dist)), method = "spearman"),
    component_relation_spearman = cor(input_relation[relation_index], code_relation[relation_index], method = "spearman")
  )
}

pca_stability <- bind_rows(pca_rows)
loading_stability <- bind_rows(loading_rows)
som_fidelity <- bind_rows(som_rows)
summary <- data.frame(
  metric = c(
    "PCA Procrustes correlation", "PCA pairwise-distance correlation",
    "PCA two-axis variance", "SOM input-to-map distance correlation",
    "SOM component-relation correlation", "SOM 10-neighbour preservation"
  ),
  median = c(
    median(pca_stability$procrustes_correlation),
    median(pca_stability$pairwise_distance_spearman),
    median(pca_stability$variance_pc1_pc2),
    median(som_fidelity$input_to_map_distance_spearman),
    median(som_fidelity$component_relation_spearman),
    median(ensemble$metrics$neighbour_preservation_10)
  ),
  lower_95 = c(
    quantile(pca_stability$procrustes_correlation, 0.025),
    quantile(pca_stability$pairwise_distance_spearman, 0.025),
    quantile(pca_stability$variance_pc1_pc2, 0.025),
    quantile(som_fidelity$input_to_map_distance_spearman, 0.025),
    quantile(som_fidelity$component_relation_spearman, 0.025),
    quantile(ensemble$metrics$neighbour_preservation_10, 0.025)
  ),
  upper_95 = c(
    quantile(pca_stability$procrustes_correlation, 0.975),
    quantile(pca_stability$pairwise_distance_spearman, 0.975),
    quantile(pca_stability$variance_pc1_pc2, 0.975),
    quantile(som_fidelity$input_to_map_distance_spearman, 0.975),
    quantile(som_fidelity$component_relation_spearman, 0.975),
    quantile(ensemble$metrics$neighbour_preservation_10, 0.975)
  )
)

gradient <- list(
  full_pca = full_pca,
  pca_stability = pca_stability,
  loading_stability = loading_stability,
  som_fidelity = som_fidelity,
  summary = summary
)
save_rds_safe(gradient, file.path(MODEL_DIR, "continuous_gradient_robustness.rds"))
write_csv_safe(pca_stability, file.path(TABLE_DIR, "TableS_pca_geometry_resampling.csv"))
write_csv_safe(loading_stability, file.path(TABLE_DIR, "TableS_pca_loading_stability.csv"))
write_csv_safe(som_fidelity, file.path(TABLE_DIR, "TableS_som_gradient_fidelity.csv"))
write_csv_safe(summary, file.path(TABLE_DIR, "TableS_continuous_gradient_summary.csv"))

cat("Continuous-gradient robustness completed\n")
print(summary)
