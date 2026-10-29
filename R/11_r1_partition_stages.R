args <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "00_config.R"))

# Attribution of partition instability to workflow stages (first-round revision).
#
# The SOM workflow has three stages: map training, Ward.D2 clustering of the
# codebook and consensus aggregation. The primary partitions clustered every
# prototype with equal weight, empty nodes included. Here the second stage is
# replaced, on the same trained maps, by Ward.D2 on occupied prototypes weighted
# by their training hits, so that any change is attributable to that stage. The
# consensus stage is checked by comparing consensus partitions built from two
# different perturbation designs. The Ward.D2 baseline is re-evaluated with
# training events keeping their tree memberships, because nearest-centroid
# allocation can move them; K-means allocation is checked in the same way.

assert_packages(c("kohonen", "mclust"))
prepared <- readRDS(file.path(RESULT_DIR, "prepared_data.rds"))
ensemble <- readRDS(file.path(MODEL_DIR, "som_ensemble_250.rds"))
resampling <- readRDS(file.path(MODEL_DIR, "cross_model_resampling_200.rds"))
resample_soms <- readRDS(file.path(MODEL_DIR, "r1_resample_soms_200.rds"))
data <- prepared$cohorts$core9
n <- nrow(data)

# Ward.D2 started from occupied prototypes treated as clusters of their hit
# counts. With d(A, B) = sqrt(2 nA nB / (nA + nB)) |cA - cB| as the starting
# dissimilarity, the Lance-Williams updates reproduce Ward.D2 on the quantised
# events (each event replaced by its prototype).
weighted_ward_tree <- function(codes, hits) {
  occupied <- which(hits > 0)
  w <- hits[occupied]
  centres <- codes[occupied, , drop = FALSE]
  start <- as.matrix(dist(centres)) * sqrt(2 * outer(w, w) / outer(w, w, "+"))
  list(tree = hclust(as.dist(start), method = "ward.D2", members = w), occupied = occupied)
}

# Check the construction once against Ward.D2 on the expanded quantised data.
check_model <- ensemble$models[[which(ensemble$metrics$grid == "8x6")[1]]]
check_hits <- tabulate(check_model$unit.classif, nbins = nrow(check_model$codes[[1]]))
check_tree <- weighted_ward_tree(check_model$codes[[1]], check_hits)
expanded_units <- rep(check_tree$occupied, check_hits[check_tree$occupied])
expanded_tree <- hclust(dist(check_model$codes[[1]][expanded_units, , drop = FALSE]), method = "ward.D2")
for (k in K_VALUES) {
  a <- cutree(check_tree$tree, k = k)
  b <- cutree(expanded_tree, k = k)[match(check_tree$occupied, expanded_units)]
  if (mclust::adjustedRandIndex(a, b) < 1 - 1e-12) stop("Weighted prototype Ward.D2 does not match expanded-data Ward.D2 at k = ", k)
}

weighted_som_labels <- function(codes, hits, k, x_target = NULL, target_bmu = NULL) {
  fit <- weighted_ward_tree(codes, hits)
  groups <- cutree(fit$tree, k = k)
  if (!is.null(target_bmu)) {
    unit_group <- integer(nrow(codes)); unit_group[fit$occupied] <- groups
    if (any(unit_group[target_bmu] == 0L)) stop("Target event mapped to an empty unit")
    return(unit_group[target_bmu])
  }
  groups[nearest_centroid(x_target, codes[fit$occupied, , drop = FALSE])]
}

stability_profile <- function(assignments, k) {
  cons <- consensus_partition(assignments, k)
  reference <- cons$labels
  fits <- seq_len(ncol(assignments))
  jaccard <- sapply(fits, function(j) clusterwise_jaccard(reference, assignments[, j])$jaccard)
  jaccard <- matrix(jaccard, nrow = length(unique(reference)))
  ari <- vapply(fits, function(j) mclust::adjustedRandIndex(reference, assignments[, j]), numeric(1))
  smallest <- vapply(fits, function(j) min(table(assignments[, j])), numeric(1))
  class_median <- as.character(stability_class(apply(jaccard, 1, median)))
  class_mean <- as.character(stability_class(rowMeans(jaccard)))
  count_classes <- function(x) sprintf("%d / %d / %d", sum(x == "stable"), sum(x == "uncertain"), sum(x == "fragile"))
  list(
    consensus = reference,
    row = data.frame(
      k = k,
      fits = length(fits),
      median_ari = median(ari),
      lower_ari = unname(quantile(ari, 0.025)),
      upper_ari = unname(quantile(ari, 0.975)),
      median_mean_jaccard = median(colMeans(jaccard)),
      classes_median_jaccard = count_classes(class_median),
      classes_mean_jaccard = count_classes(class_mean),
      fragile_clusters = sum(class_median == "fragile"),
      consensus_sizes = paste(sort(as.integer(table(reference)), decreasing = TRUE), collapse = "/"),
      median_smallest_cluster = median(smallest),
      share_fits_smallest_le5 = mean(smallest <= 5)
    )
  )
}

# SOM second-level clustering: unweighted (primary) versus hit-weighted ---------
som_rows <- list(); som_consensus <- list(); counter <- 0L
for (k in K_VALUES) {
  key <- paste0("k", k)
  grid_seed_weighted <- sapply(ensemble$models, function(model) {
    codes <- model$codes[[1]]
    weighted_som_labels(codes, tabulate(model$unit.classif, nbins = nrow(codes)), k, target_bmu = model$unit.classif)
  })
  resample_weighted <- matrix(NA_integer_, n, length(resample_soms))
  for (b in seq_along(resample_soms)) {
    fit <- resample_soms[[b]]
    train <- data[fit$train_idx, , drop = FALSE]
    x_all <- transform_apply(data, transform_fit(train))
    hits <- tabulate(fit$train_bmu, nbins = nrow(fit$codes))
    resample_weighted[, b] <- weighted_som_labels(fit$codes, hits, k, x_target = x_all)
  }
  variants <- list(
    list(design = "Grid-seed ensemble (250 full-data maps)", variant = "Unweighted, all prototypes (primary)", assignments = ensemble$assignments[[key]]),
    list(design = "Grid-seed ensemble (250 full-data maps)", variant = "Hit-weighted, occupied prototypes", assignments = grid_seed_weighted),
    list(design = "Grouped resampling (200 maps)", variant = "Unweighted, all prototypes (primary)", assignments = resampling$assignments$SOM[[key]]),
    list(design = "Grouped resampling (200 maps)", variant = "Hit-weighted, occupied prototypes", assignments = resample_weighted)
  )
  for (v in variants) {
    profile <- stability_profile(v$assignments, k)
    counter <- counter + 1L
    som_consensus[[paste(v$design, v$variant, key)]] <- profile$consensus
    som_rows[[counter]] <- profile$row %>%
      mutate(
        design = v$design, variant = v$variant,
        ari_to_primary_consensus = mclust::adjustedRandIndex(profile$consensus, resampling$consensus$SOM[[key]]$labels),
        .before = 1
      )
  }
  message("SOM second-level variants: k = ", k)
}
som_second_level <- bind_rows(som_rows)

consensus_level <- bind_rows(lapply(K_VALUES, function(k) {
  key <- paste0("k", k)
  get <- function(design, variant) som_consensus[[paste(design, variant, key)]]
  data.frame(
    k = k,
    unweighted_grid_seed_vs_resampling = mclust::adjustedRandIndex(
      get("Grid-seed ensemble (250 full-data maps)", "Unweighted, all prototypes (primary)"),
      get("Grouped resampling (200 maps)", "Unweighted, all prototypes (primary)")),
    weighted_grid_seed_vs_resampling = mclust::adjustedRandIndex(
      get("Grid-seed ensemble (250 full-data maps)", "Hit-weighted, occupied prototypes"),
      get("Grouped resampling (200 maps)", "Hit-weighted, occupied prototypes")),
    weighted_vs_unweighted_resampling = mclust::adjustedRandIndex(
      get("Grouped resampling (200 maps)", "Hit-weighted, occupied prototypes"),
      get("Grouped resampling (200 maps)", "Unweighted, all prototypes (primary)"))
  )
}))

# Ward.D2 and K-means allocation of training events ----------------------------
ward_tree_assign <- setNames(lapply(K_VALUES, function(k) matrix(NA_integer_, n, N_RESAMPLE)), paste0("k", K_VALUES))
allocation_rows <- list(); counter <- 0L
for (b in seq_len(N_RESAMPLE)) {
  train_idx <- stratified_block_sample(data, RESAMPLE_FRACTION, RESAMPLE_SEEDS[b])
  train <- data[train_idx, , drop = FALSE]
  transform <- transform_fit(train)
  x_train <- transform_apply(train, transform)
  x_all <- transform_apply(data, transform)
  tree <- hclust(dist(x_train), method = "ward.D2")
  for (k in K_VALUES) {
    tree_labels <- cutree(tree, k = k)
    centres <- do.call(rbind, lapply(sort(unique(tree_labels)), function(g) colMeans(x_train[tree_labels == g, , drop = FALSE])))
    nearest_train <- nearest_centroid(x_train, centres)
    labels <- nearest_centroid(x_all, centres)
    labels[train_idx] <- tree_labels
    ward_tree_assign[[paste0("k", k)]][, b] <- labels

    set.seed(RESAMPLE_SEEDS[b] + k)
    km <- kmeans(x_train, centers = k, nstart = 50, iter.max = 200)
    counter <- counter + 1L
    allocation_rows[[counter]] <- data.frame(
      resample = b, k = k,
      ward_tree_vs_nearest_ari = mclust::adjustedRandIndex(tree_labels, nearest_train),
      ward_tree_vs_nearest_identical = mean(tree_labels == nearest_train),
      kmeans_fit_vs_nearest_identical = mean(km$cluster == nearest_centroid(x_train, km$centers))
    )
  }
  if (b %% 50L == 0L) message("Ward and K-means allocation checks: ", b, "/", N_RESAMPLE)
}
allocation <- bind_rows(allocation_rows)

ward_rows <- lapply(K_VALUES, function(k) {
  key <- paste0("k", k)
  profile <- stability_profile(ward_tree_assign[[key]], k)
  alloc <- allocation[allocation$k == k, ]
  profile$row %>% mutate(
    ward_tree_vs_nearest_median_ari = median(alloc$ward_tree_vs_nearest_ari),
    ward_tree_vs_nearest_min_ari = min(alloc$ward_tree_vs_nearest_ari),
    ward_tree_vs_nearest_median_identical = median(alloc$ward_tree_vs_nearest_identical),
    ward_tree_vs_nearest_min_identical = min(alloc$ward_tree_vs_nearest_identical),
    kmeans_share_fits_identical = mean(alloc$kmeans_fit_vs_nearest_identical == 1),
    kmeans_min_identical = min(alloc$kmeans_fit_vs_nearest_identical),
    ari_tree_consensus_vs_primary_ward = mclust::adjustedRandIndex(profile$consensus, resampling$consensus[["Ward.D2"]][[key]]$labels),
    ari_tree_consensus_vs_kmeans = mclust::adjustedRandIndex(profile$consensus, resampling$consensus[["K-means"]][[key]]$labels),
    ari_primary_ward_vs_kmeans = mclust::adjustedRandIndex(resampling$consensus[["Ward.D2"]][[key]]$labels, resampling$consensus[["K-means"]][[key]]$labels),
    .before = 1
  )
})
ward_allocation <- bind_rows(ward_rows)

stages <- list(
  som_second_level = som_second_level,
  consensus_level = consensus_level,
  som_consensus = som_consensus,
  ward_allocation = ward_allocation,
  ward_tree_assignments = ward_tree_assign,
  allocation_by_resample = allocation
)
save_rds_safe(stages, file.path(MODEL_DIR, "r1_partition_stages.rds"))
write_csv_safe(som_second_level, file.path(TABLE_DIR, "TableS_r1_som_second_level.csv"))
write_csv_safe(consensus_level, file.path(TABLE_DIR, "TableS_r1_consensus_level.csv"))
write_csv_safe(ward_allocation, file.path(TABLE_DIR, "TableS_r1_ward_kmeans_allocation.csv"))

cat("Partition-stage checks completed\n")
print(som_second_level %>% filter(k %in% 3:4) %>% select(design, variant, k, median_ari, classes_median_jaccard, consensus_sizes, ari_to_primary_consensus))
print(consensus_level)
print(ward_allocation %>% select(k, median_ari, classes_median_jaccard, ward_tree_vs_nearest_median_identical, ari_tree_consensus_vs_primary_ward, ari_tree_consensus_vs_kmeans))
