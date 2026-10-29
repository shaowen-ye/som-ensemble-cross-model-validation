args <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "00_config.R"))

# Definitions and summaries of the stability metrics (first-round revision).
#
# 1. Clusterwise Jaccard classes from the mean over resamples (Hennig's guidelines
#    refer to mean bootstrap Jaccard) beside the median classes used in the
#    submitted analysis, and median classes from training events only.
# 2. Sample support, margin and entropy with the event's own coassociation (1)
#    excluded from its cluster mean; singletons have no support under this
#    definition and are flagged.
# 3. Median, interquartile range and mean of normalised entropy for every
#    method and k, so that the figures can be checked against one table.
# 4. The training membership of every grouped resample, with the pseudonymised
#    event and block identifiers of the public data file.

assert_packages(c("mclust"))
prepared <- readRDS(file.path(RESULT_DIR, "prepared_data.rds"))
resampling <- readRDS(file.path(MODEL_DIR, "cross_model_resampling_200.rds"))
data <- prepared$cohorts$core9
n <- nrow(data)
methods <- c("SOM", "K-means", "Ward.D2", "GMM")

# 1. Mean versus median clusterwise Jaccard ----------------------------------
jaccard_rows <- list(); class_rows <- list(); counter <- 0L
for (method in methods) for (k in K_VALUES) {
  key <- paste0("k", k)
  reference <- resampling$consensus[[method]][[key]]$labels
  values <- sapply(seq_len(ncol(resampling$assignments[[method]][[key]])), function(b) {
    clusterwise_jaccard(reference, resampling$assignments[[method]][[key]][, b])$jaccard
  })
  values <- matrix(values, nrow = length(unique(reference)))
  clusters <- sort(unique(reference))
  med <- apply(values, 1, median); avg <- rowMeans(values)
  counter <- counter + 1L
  jaccard_rows[[counter]] <- data.frame(
    method = method, k = k, cluster = clusters,
    cluster_size = as.integer(table(reference)[as.character(clusters)]),
    median_jaccard = med, mean_jaccard = avg,
    class_median = as.character(stability_class(med)),
    class_mean = as.character(stability_class(avg))
  )
  count_classes <- function(x) sprintf("%d / %d / %d", sum(x == "stable"), sum(x == "uncertain"), sum(x == "fragile"))
  class_rows[[counter]] <- data.frame(
    method = method, k = k,
    classes_median = count_classes(as.character(stability_class(med))),
    classes_mean = count_classes(as.character(stability_class(avg))),
    clusters_changed = sum(as.character(stability_class(med)) != as.character(stability_class(avg))),
    fragile_median = sum(med < 0.60), fragile_mean = sum(avg < 0.60)
  )
}
jaccard_by_cluster <- bind_rows(jaccard_rows)
jaccard_classes <- bind_rows(class_rows)

# 1b. Jaccard over each resample's training events only -------------------------
# Hennig's bootstrap Jaccard is computed on the resampled points; the primary
# values use all 390 events, about a fifth of which are predicted outside the
# training set. Clusters absent from a training set are skipped for that resample.
train_sets <- lapply(seq_len(N_RESAMPLE), function(b) stratified_block_sample(data, RESAMPLE_FRACTION, RESAMPLE_SEEDS[b]))
training_rows <- list(); counter <- 0L
for (method in methods) for (k in K_VALUES) {
  key <- paste0("k", k)
  reference <- resampling$consensus[[method]][[key]]$labels
  clusters <- sort(unique(reference))
  values <- sapply(seq_len(N_RESAMPLE), function(b) {
    idx <- train_sets[[b]]
    j <- clusterwise_jaccard(reference[idx], resampling$assignments[[method]][[key]][idx, b])
    out <- rep(NA_real_, length(clusters))
    out[match(j$cluster, clusters)] <- j$jaccard
    out
  })
  values <- matrix(values, nrow = length(clusters))
  med <- apply(values, 1, median, na.rm = TRUE)
  primary <- jaccard_by_cluster$class_median[jaccard_by_cluster$method == method & jaccard_by_cluster$k == k]
  counter <- counter + 1L
  training_rows[[counter]] <- data.frame(
    method = method, k = k,
    classes_training_only = count_classes(as.character(stability_class(med))),
    clusters_changed_training_only = sum(as.character(stability_class(med)) != primary),
    fragile_training_only = sum(med < 0.60)
  )
}
jaccard_classes <- jaccard_classes %>% left_join(bind_rows(training_rows), by = c("method", "k"))

# 2. Sample metrics with self-pairs excluded ----------------------------------
sample_metrics_excluding_self <- function(labels, coassoc) {
  groups <- sort(unique(labels))
  rows <- lapply(seq_along(labels), function(i) {
    affinities <- vapply(groups, function(g) {
      members <- which(labels == g)
      members <- members[members != i]
      if (!length(members)) NA_real_ else mean(coassoc[i, members])
    }, numeric(1))
    own <- affinities[groups == labels[i]]
    singleton <- is.na(own)
    probabilities <- affinities / sum(affinities, na.rm = TRUE)
    data.frame(
      row_index = i,
      singleton = singleton,
      support_excluding_self = if (singleton) NA_real_ else own,
      margin_excluding_self = if (singleton) NA_real_ else own - max(affinities[groups != labels[i]], na.rm = TRUE),
      entropy_excluding_self = if (singleton) NA_real_ else entropy(probabilities[!is.na(probabilities)]) / log(length(groups))
    )
  })
  bind_rows(rows)
}

sample_rows <- list(); counter <- 0L
for (method in methods) for (k in K_VALUES) {
  key <- paste0("k", k)
  cons <- resampling$consensus[[method]][[key]]
  original <- resampling$sample_metrics %>% filter(method == !!method, k == !!k) %>% arrange(row_index)
  excluded <- sample_metrics_excluding_self(cons$labels, cons$coassociation)
  counter <- counter + 1L
  sample_rows[[counter]] <- original %>%
    select(method, k, row_index, consensus_cluster, support, margin, entropy) %>%
    left_join(excluded, by = "row_index")
}
sample_metrics <- bind_rows(sample_rows)

entropy_summary <- sample_metrics %>%
  group_by(method, k) %>%
  summarise(
    events = n(),
    singleton_events = sum(singleton),
    median_entropy = median(entropy),
    q25_entropy = unname(quantile(entropy, 0.25)),
    q75_entropy = unname(quantile(entropy, 0.75)),
    mean_entropy = mean(entropy),
    median_entropy_excluding_self = median(entropy_excluding_self, na.rm = TRUE),
    median_support = median(support),
    median_support_excluding_self = median(support_excluding_self, na.rm = TRUE),
    median_margin = median(margin),
    median_margin_excluding_self = median(margin_excluding_self, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(method = factor(method, levels = methods)) %>%
  arrange(method, k) %>%
  mutate(method = as.character(method))

entropy_rank <- entropy_summary %>%
  group_by(k) %>%
  summarise(order_by_median_entropy = paste(method[order(-median_entropy)], collapse = " > "), .groups = "drop")

# 4. Training membership of the grouped resamples -----------------------------
membership <- matrix(0L, n, N_RESAMPLE)
for (b in seq_len(N_RESAMPLE)) {
  membership[stratified_block_sample(data, RESAMPLE_FRACTION, RESAMPLE_SEEDS[b]), b] <- 1L
}
if (!identical(as.integer(colSums(membership)), as.integer(resampling$resample_manifest$n_train))) {
  stop("Exported training membership does not match the resampling manifest")
}
colnames(membership) <- sprintf("R%03d", seq_len(N_RESAMPLE))
# Event and block identifiers in the public data are already pseudonymised.
index_table <- data.frame(event_id = data$record_id, block_id_public = data$block_id,
                          campaign = as.character(data$campaign))
index_table <- bind_cols(index_table %>% select(event_id, block_id_public, campaign), as.data.frame(membership))
resample_seeds <- data.frame(resample = sprintf("R%03d", seq_len(N_RESAMPLE)), seed = RESAMPLE_SEEDS,
                             som_grid = resampling$resample_manifest$som_grid, n_train = resampling$resample_manifest$n_train)

checks <- list(
  jaccard_by_cluster = jaccard_by_cluster,
  jaccard_classes = jaccard_classes,
  sample_metrics = sample_metrics,
  entropy_summary = entropy_summary,
  entropy_rank = entropy_rank
)
save_rds_safe(checks, file.path(MODEL_DIR, "r1_metric_checks.rds"))
write_csv_safe(jaccard_by_cluster, file.path(TABLE_DIR, "TableS_r1_jaccard_mean_vs_median_by_cluster.csv"))
write_csv_safe(jaccard_classes, file.path(TABLE_DIR, "TableS_r1_jaccard_mean_vs_median_classes.csv"))
write_csv_safe(entropy_summary, file.path(TABLE_DIR, "TableS_r1_sample_metric_summary.csv"))
write_csv_safe(entropy_rank, file.path(TABLE_DIR, "TableS_r1_entropy_order_by_k.csv"))
write_csv_safe(index_table, file.path(TABLE_DIR, "TableS_r1_resampling_training_membership_200.csv"))
write_csv_safe(resample_seeds, file.path(TABLE_DIR, "TableS_r1_resampling_seeds_200.csv"))

cat("Metric checks completed\n")
print(jaccard_classes %>% filter(clusters_changed > 0))
print(entropy_summary %>% filter(k == 3L) %>% select(method, median_entropy, q25_entropy, q75_entropy, mean_entropy, median_entropy_excluding_self, singleton_events))
print(entropy_rank)
