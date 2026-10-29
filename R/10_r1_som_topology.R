args <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "00_config.R"))

# Reproducibility of the SOM map layout (first-round revision).
#
# Component-relation correlation is unchanged when prototypes are rearranged on
# the grid, and PCA agreement says nothing about the map itself. This script
# compares maps directly. For every pair of SOM fits it correlates the grid
# distances between the best matching units (BMUs) of all event pairs, and it
# measures how far each event keeps the same map neighbours (events whose BMU is
# the same or an adjacent unit). A random relabelling of events gives the chance
# level of the neighbour overlap. The 200 grouped-resample SOMs of
# 04_resampling_validation.R are refitted with the same seeds and grids; their
# partitions must match the stored ones exactly.
#
# Grid distances between BMUs are discrete and heavily tied, so the continuous
# PCA plane is not a commensurate reference. As one, the aligned PCA scores of
# each resample are quantised onto the hexagonal grid of the matching SOM
# (marginal ranks scaled to the grid, nearest unit) and compared in the same
# way. Component-relation correlation is also computed for the resample maps.

assert_packages(c("kohonen", "mclust"))
prepared <- readRDS(file.path(RESULT_DIR, "prepared_data.rds"))
ensemble <- readRDS(file.path(MODEL_DIR, "som_ensemble_250.rds"))
resampling <- readRDS(file.path(MODEL_DIR, "cross_model_resampling_200.rds"))
data <- prepared$cohorts$core9
n <- nrow(data)
upper <- upper.tri(matrix(FALSE, n, n))

map_distance_matrix <- function(grid, bmu) {
  kohonen::unit.distances(grid)[bmu, bmu, drop = FALSE]
}

neighbour_matrix <- function(map_distance) {
  out <- map_distance <= 1 + 1e-8
  diag(out) <- FALSE
  out
}

neighbour_overlap <- function(a, b) {
  union <- rowSums(a | b)
  keep <- union > 0
  mean(rowSums(a & b)[keep] / union[keep])
}

align_to_reference <- function(scores, reference) {
  s <- svd(crossprod(scores, reference))
  scores %*% (s$u %*% t(s$v))
}

quantise_to_grid <- function(scores, grid) {
  pts <- grid$pts
  u <- apply(scores, 2, function(z) (rank(z, ties.method = "average") - 0.5) / length(z))
  target <- cbind(min(pts[, 1]) + u[, 1] * diff(range(pts[, 1])),
                  min(pts[, 2]) + u[, 2] * diff(range(pts[, 2])))
  nearest_units(target, pts)$unit
}

between_fit_concordance <- function(distance_vectors) {
  ranks <- apply(distance_vectors, 2, rank)
  stats::cor(ranks)
}

pair_table <- function(concordance, labels, design) {
  index <- which(upper.tri(concordance), arr.ind = TRUE)
  data.frame(
    design = design,
    fit_1 = labels$fit[index[, 1]],
    fit_2 = labels$fit[index[, 2]],
    grid_1 = labels$grid[index[, 1]],
    grid_2 = labels$grid[index[, 2]],
    map_distance_spearman = concordance[index]
  )
}

overlap_rows <- function(neighbours, labels, design, seed) {
  set.seed(seed)
  rows <- list(); counter <- 0L
  for (grid in unique(labels$grid)) {
    members <- which(labels$grid == grid)
    if (length(members) < 2L) next
    pairs <- combn(members, 2)
    for (j in seq_len(ncol(pairs))) {
      a <- neighbours[[pairs[1, j]]]
      b <- neighbours[[pairs[2, j]]]
      permutation <- sample.int(n)
      counter <- counter + 1L
      rows[[counter]] <- data.frame(
        design = design,
        grid = grid,
        fit_1 = labels$fit[pairs[1, j]],
        fit_2 = labels$fit[pairs[2, j]],
        neighbour_overlap = neighbour_overlap(a, b),
        chance_overlap = neighbour_overlap(a, b[permutation, permutation])
      )
    }
  }
  bind_rows(rows)
}

# Full-data ensemble: 250 maps (five grids x 50 seeds) -----------------------
full_labels <- data.frame(fit = ensemble$metrics$model, grid = ensemble$metrics$grid)
full_vectors <- matrix(NA_real_, sum(upper), length(ensemble$models))
full_neighbours <- vector("list", length(ensemble$models))
for (i in seq_along(ensemble$models)) {
  model <- ensemble$models[[i]]
  md <- map_distance_matrix(model$grid, model$unit.classif)
  full_vectors[, i] <- md[upper]
  full_neighbours[[i]] <- neighbour_matrix(md)
}
full_concordance <- between_fit_concordance(full_vectors)
full_pairs <- pair_table(full_concordance, full_labels, "Full-data SOM ensemble")
full_overlap <- overlap_rows(full_neighbours, full_labels, "Full-data SOM ensemble", 97001L)
rm(full_vectors)
message("Full-data map comparison completed")

# Grouped-resample SOMs: refit with the seeds and grids of 04 ------------------
resample_soms <- vector("list", N_RESAMPLE)
resample_vectors <- matrix(NA_real_, sum(upper), N_RESAMPLE)
resample_neighbours <- vector("list", N_RESAMPLE)
pca_vectors <- matrix(NA_real_, sum(upper), N_RESAMPLE)
quantised_vectors <- matrix(NA_real_, sum(upper), N_RESAMPLE)
resample_relation <- numeric(N_RESAMPLE)
relation_index <- upper.tri(diag(length(VARS)))
x_full <- transform_apply(data, transform_fit(data))
full_scores <- x_full %*% prcomp(x_full, center = FALSE, scale. = FALSE)$rotation[, 1:2, drop = FALSE]
mismatch_rows <- list()
for (b in seq_len(N_RESAMPLE)) {
  train_idx <- stratified_block_sample(data, RESAMPLE_FRACTION, RESAMPLE_SEEDS[b])
  train <- data[train_idx, , drop = FALSE]
  transform <- transform_fit(train)
  x_train <- transform_apply(train, transform)
  x_all <- transform_apply(data, transform)
  g <- 1L + ((b - 1L) %% nrow(GRIDS))
  set.seed(RESAMPLE_SEEDS[b])
  som_fit <- kohonen::som(
    x_train,
    grid = kohonen::somgrid(GRIDS$xdim[g], GRIDS$ydim[g], topo = "hexagonal"),
    rlen = SOM_RLEN,
    alpha = c(0.05, 0.01),
    keep.data = TRUE,
    mode = "online"
  )
  mismatch_rows[[b]] <- data.frame(
    resample = b,
    k = K_VALUES,
    mismatched_labels = vapply(K_VALUES, function(k) {
      sum(partition_som(som_fit, k, x_all) != resampling$assignments$SOM[[paste0("k", k)]][, b])
    }, numeric(1))
  )
  all_bmu <- nearest_units(x_all, som_fit$codes[[1]])$unit
  md <- map_distance_matrix(som_fit$grid, all_bmu)
  resample_vectors[, b] <- md[upper]
  resample_neighbours[[b]] <- neighbour_matrix(md)
  resample_soms[[b]] <- list(
    resample = b,
    seed = RESAMPLE_SEEDS[b],
    xdim = GRIDS$xdim[g],
    ydim = GRIDS$ydim[g],
    train_idx = train_idx,
    codes = som_fit$codes[[1]],
    train_bmu = som_fit$unit.classif,
    all_bmu = all_bmu
  )

  resample_relation[b] <- cor(cor(x_train)[relation_index], cor(som_fit$codes[[1]])[relation_index], method = "spearman")

  pca <- prcomp(x_train, center = FALSE, scale. = FALSE)
  scores <- x_all %*% pca$rotation[, 1:2, drop = FALSE]
  pca_vectors[, b] <- as.numeric(dist(scores))
  cells <- quantise_to_grid(align_to_reference(scores, full_scores), som_fit$grid)
  quantised_vectors[, b] <- map_distance_matrix(som_fit$grid, cells)[upper]
  if (b %% 20L == 0L) message("Grouped-resample maps: ", b, "/", N_RESAMPLE)
}
mismatches <- bind_rows(mismatch_rows)
if (sum(mismatches$mismatched_labels) != 0) {
  stop("Refitted resample SOMs do not reproduce the stored partitions")
}

resample_labels <- data.frame(
  fit = sprintf("resample_%03d", seq_len(N_RESAMPLE)),
  grid = vapply(resample_soms, function(z) sprintf("%dx%d", z$xdim, z$ydim), character(1))
)
resample_concordance <- between_fit_concordance(resample_vectors)
resample_pairs <- pair_table(resample_concordance, resample_labels, "Grouped-resample SOMs")
resample_overlap <- overlap_rows(resample_neighbours, resample_labels, "Grouped-resample SOMs", 97002L)
rm(resample_vectors)

pca_concordance <- between_fit_concordance(pca_vectors)
pca_pairs <- data.frame(
  design = "Grouped-resample PCA (two axes)",
  map_distance_spearman = pca_concordance[upper.tri(pca_concordance)]
)
rm(pca_vectors)

quantised_concordance <- between_fit_concordance(quantised_vectors)
quantised_pairs <- pair_table(quantised_concordance, resample_labels, "Grouped-resample PCA quantised to the SOM grid")
rm(quantised_vectors)

# Summaries ----------------------------------------------------------------
summarise_metric <- function(x, design, comparison, metric, grid = "All grids") {
  data.frame(
    design = design, comparison = comparison, grid = grid, metric = metric,
    pairs = length(x),
    median = median(x), lower_95 = unname(quantile(x, 0.025)), upper_95 = unname(quantile(x, 0.975))
  )
}

summary_rows <- list(
  summarise_metric(full_pairs$map_distance_spearman, "Full-data SOM ensemble", "All pairs", "Map-distance concordance"),
  summarise_metric(full_pairs$map_distance_spearman[full_pairs$grid_1 == full_pairs$grid_2], "Full-data SOM ensemble", "Same-grid pairs", "Map-distance concordance"),
  summarise_metric(full_overlap$neighbour_overlap, "Full-data SOM ensemble", "Same-grid pairs", "Map-neighbour overlap"),
  summarise_metric(full_overlap$chance_overlap, "Full-data SOM ensemble", "Same-grid pairs", "Map-neighbour overlap, chance"),
  summarise_metric(resample_pairs$map_distance_spearman, "Grouped-resample SOMs", "All pairs", "Map-distance concordance"),
  summarise_metric(resample_pairs$map_distance_spearman[resample_pairs$grid_1 == resample_pairs$grid_2], "Grouped-resample SOMs", "Same-grid pairs", "Map-distance concordance"),
  summarise_metric(resample_overlap$neighbour_overlap, "Grouped-resample SOMs", "Same-grid pairs", "Map-neighbour overlap"),
  summarise_metric(resample_overlap$chance_overlap, "Grouped-resample SOMs", "Same-grid pairs", "Map-neighbour overlap, chance"),
  summarise_metric(resample_relation, "Grouped-resample SOMs", "Single maps", "Component-relation correlation"),
  summarise_metric(pca_pairs$map_distance_spearman, "Grouped-resample PCA (two axes)", "All pairs", "Score-distance concordance"),
  summarise_metric(quantised_pairs$map_distance_spearman, "Grouped-resample PCA quantised to the SOM grid", "All pairs", "Map-distance concordance"),
  summarise_metric(quantised_pairs$map_distance_spearman[quantised_pairs$grid_1 == quantised_pairs$grid_2], "Grouped-resample PCA quantised to the SOM grid", "Same-grid pairs", "Map-distance concordance")
)
for (grid in levels(factor(full_labels$grid))) {
  same <- full_pairs$grid_1 == grid & full_pairs$grid_2 == grid
  summary_rows[[length(summary_rows) + 1L]] <- summarise_metric(full_pairs$map_distance_spearman[same], "Full-data SOM ensemble", "Same-grid pairs", "Map-distance concordance", grid)
  summary_rows[[length(summary_rows) + 1L]] <- summarise_metric(full_overlap$neighbour_overlap[full_overlap$grid == grid], "Full-data SOM ensemble", "Same-grid pairs", "Map-neighbour overlap", grid)
}
topology_summary <- bind_rows(summary_rows)

topology <- list(
  summary = topology_summary,
  full_pairs = full_pairs,
  full_overlap = full_overlap,
  resample_pairs = resample_pairs,
  resample_overlap = resample_overlap,
  pca_pairs = pca_pairs,
  quantised_pairs = quantised_pairs,
  resample_component_relation = resample_relation,
  refit_mismatches = mismatches
)
save_rds_safe(topology, file.path(MODEL_DIR, "r1_som_topology.rds"))
save_rds_safe(resample_soms, file.path(MODEL_DIR, "r1_resample_soms_200.rds"))
write_csv_safe(topology_summary, file.path(TABLE_DIR, "TableS_r1_som_map_reproducibility.csv"))
write_csv_safe(mismatches %>% group_by(k) %>% summarise(mismatched_labels = sum(mismatched_labels), .groups = "drop"),
               file.path(TABLE_DIR, "TableS_r1_resample_som_refit_check.csv"))

cat("Refitted resample SOMs reproduce stored partitions: all", N_RESAMPLE, "maps, k = 2-8\n")
print(topology_summary %>% filter(grid == "All grids"))
