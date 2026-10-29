args <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "00_config.R"))

# Calibration of survey transfer against matched-size random omission
# (first-round revision).
#
# Holding out one survey changes both the environmental composition of the
# training data and its size (16 to 186 omitted events). Each leave-one-survey-out
# (LOSO) fit of 05_campaign_transfer.R is therefore compared with fits in which
# the same number of events is removed at random: site-within-survey blocks are
# drawn from all surveys until the omitted count is closest to the held-out
# survey size. Fitting, prediction, the concordance reference and the
# centroid-distance ratio follow 05 exactly. SOM agreement is reported both as
# the median over its ten single maps and for the consensus of those ten maps.
# A single survey can be more homogeneous in its reference labels than a random
# omission of the same size, and ARI is sensitive to small reference clusters,
# so AMI is reported beside ARI and the number of reference clusters and the
# share of the largest one are recorded for every method and test set.

assert_packages(c("kohonen", "mclust", "parallel"))
prepared <- readRDS(file.path(RESULT_DIR, "prepared_data.rds"))
resampling <- readRDS(file.path(MODEL_DIR, "cross_model_resampling_200.rds"))
transfer <- readRDS(file.path(MODEL_DIR, "campaign_transfer.rds"))
data <- prepared$cohorts$core9
methods <- c("SOM", "K-means", "Ward.D2", "GMM")
cores <- as.integer(Sys.getenv("EI_CORES", unset = as.character(max(1L, min(4L, parallel::detectCores() - 1L)))))

centroid_distance_ratio <- function(x_train, x_test, train_labels, test_labels) {
  groups <- sort(unique(train_labels))
  centers <- do.call(rbind, lapply(groups, function(g) colMeans(x_train[train_labels == g, , drop = FALSE])))
  rownames(centers) <- groups
  train_distance <- vapply(seq_len(nrow(x_train)), function(i) sqrt(sum((x_train[i, ] - centers[as.character(train_labels[i]), ])^2)), numeric(1))
  test_distance <- vapply(seq_len(nrow(x_test)), function(i) {
    label <- as.character(test_labels[i])
    if (!label %in% rownames(centers)) return(NA_real_)
    sqrt(sum((x_test[i, ] - centers[label, ])^2))
  }, numeric(1))
  median(test_distance, na.rm = TRUE) / median(train_distance, na.rm = TRUE)
}

matched_omission <- function(target_n, seed) {
  set.seed(seed)
  # Blocks are listed in order of first appearance, so the draw does not depend
  # on how block identifiers are spelled or sorted in the current locale.
  blocks <- data.frame(block_id = unique(data$block_id), stringsAsFactors = FALSE)
  blocks$events <- as.integer(table(data$block_id)[blocks$block_id])
  blocks <- blocks[sample.int(nrow(blocks)), , drop = FALSE]
  cut <- which.min(abs(cumsum(blocks$events) - target_n))
  which(data$block_id %in% blocks$block_id[seq_len(cut)])
}

held_out_fit <- function(test_idx) {
  train_idx <- setdiff(seq_len(nrow(data)), test_idx)
  train <- data[train_idx, , drop = FALSE]
  test <- data[test_idx, , drop = FALSE]
  transform <- transform_fit(train)
  x_train <- transform_apply(train, transform)
  x_test <- transform_apply(test, transform)
  x_combined <- rbind(x_train, x_test)
  n_train <- nrow(x_train)
  split_labels <- function(labels) list(train = labels[seq_len(n_train)], test = labels[n_train + seq_len(nrow(x_test))])
  score <- function(method, repeat_id, k, reference, s) {
    data.frame(method = method, repeat_id = repeat_id, k = k,
      ari = if (is.null(s)) NA_real_ else mclust::adjustedRandIndex(reference, s$test),
      ami = if (is.null(s)) NA_real_ else adjusted_mutual_information(reference, s$test),
      distance_ratio = if (is.null(s)) NA_real_ else centroid_distance_ratio(x_train, x_test, s$train, s$test),
      reference_clusters = length(unique(reference)),
      reference_largest_share = max(table(reference)) / length(reference))
  }
  rows <- list(); counter <- 0L
  som_test <- setNames(lapply(K_VALUES, function(k) matrix(NA_integer_, nrow(x_test), length(TRANSFER_SEEDS))), paste0("k", K_VALUES))

  for (repeat_id in seq_along(TRANSFER_SEEDS)) {
    seed <- TRANSFER_SEEDS[repeat_id]
    g <- 1L + ((repeat_id - 1L) %% nrow(GRIDS))
    set.seed(seed)
    som <- kohonen::som(
      x_train,
      grid = kohonen::somgrid(GRIDS$xdim[g], GRIDS$ydim[g], topo = "hexagonal"),
      rlen = SOM_RLEN, alpha = c(0.05, 0.01), keep.data = TRUE, mode = "online"
    )
    for (k in K_VALUES) {
      key <- paste0("k", k)
      s <- split_labels(partition_som(som, k, x_combined))
      som_test[[key]][, repeat_id] <- s$test
      counter <- counter + 1L
      rows[[counter]] <- score("SOM", repeat_id, k, resampling$consensus$SOM[[key]]$labels[test_idx], s)
      s <- split_labels(fit_kmeans_predict(x_train, x_combined, k, seed + k)$labels)
      counter <- counter + 1L
      rows[[counter]] <- score("K-means", repeat_id, k, resampling$consensus[["K-means"]][[key]]$labels[test_idx], s)
    }
  }
  for (k in K_VALUES) {
    key <- paste0("k", k)
    s <- split_labels(fit_ward_predict(x_train, x_combined, k)$labels)
    counter <- counter + 1L
    rows[[counter]] <- score("Ward.D2", 1L, k, resampling$consensus[["Ward.D2"]][[key]]$labels[test_idx], s)
    gmm <- fit_gmm_predict(x_train, x_combined, k)
    counter <- counter + 1L
    rows[[counter]] <- score("GMM", 1L, k, resampling$consensus$GMM[[key]]$labels[test_idx],
                             if (is.na(gmm$failure)) split_labels(gmm$labels) else NULL)
  }
  fits <- bind_rows(rows)
  consensus_rows <- lapply(K_VALUES, function(k) {
    key <- paste0("k", k)
    reference <- resampling$consensus$SOM[[key]]$labels[test_idx]
    data.frame(k = k,
      som_consensus_ari = mclust::adjustedRandIndex(reference, consensus_partition(som_test[[key]], k)$labels))
  })
  summary <- fits %>%
    group_by(method, k) %>%
    summarise(ari = median(ari), ami = median(ami), distance_ratio = median(distance_ratio),
              reference_clusters = first(reference_clusters), reference_largest_share = first(reference_largest_share),
              .groups = "drop") %>%
    left_join(bind_rows(consensus_rows), by = "k") %>%
    mutate(som_consensus_ari = ifelse(method == "SOM", som_consensus_ari, NA_real_))
  list(fits = fits, summary = summary, n_test = length(test_idx),
       surveys_in_test = length(unique(data$campaign[test_idx])))
}

# Leave-one-survey-out, refitted with stored labels ----------------------------
loso <- lapply(CAMPAIGNS, function(campaign_out) held_out_fit(which(data$campaign == campaign_out)))
names(loso) <- CAMPAIGNS
loso_fits <- bind_rows(lapply(CAMPAIGNS, function(s) loso[[s]]$fits %>% mutate(held_out_campaign = s, .before = 1)))

stored <- transfer$raw %>%
  transmute(held_out_campaign, method, repeat_id, k, stored_ari = ari_to_full_consensus,
            stored_ami = ami_to_full_consensus, stored_ratio = distance_ratio)
check <- loso_fits %>% inner_join(stored, by = c("held_out_campaign", "method", "repeat_id", "k"))
if (nrow(check) != nrow(stored) || max(abs(check$ari - check$stored_ari), na.rm = TRUE) > 1e-10 ||
    max(abs(check$ami - check$stored_ami), na.rm = TRUE) > 1e-10 ||
    max(abs(check$distance_ratio - check$stored_ratio), na.rm = TRUE) > 1e-10) {
  stop("Refitted LOSO fits do not reproduce 05_campaign_transfer.R")
}
message("LOSO refit reproduces the stored transfer results")

loso_summary <- bind_rows(lapply(CAMPAIGNS, function(s) loso[[s]]$summary %>%
  mutate(held_out_campaign = s, n_test = loso[[s]]$n_test, .before = 1)))

# Matched-size random omission -------------------------------------------------
jobs <- expand.grid(replicate = seq_len(MATCHED_OMISSION_REPLICATES), survey = seq_along(CAMPAIGNS))
run_job <- function(j) {
  survey <- CAMPAIGNS[jobs$survey[j]]
  target <- sum(data$campaign == survey)
  seed <- MATCHED_OMISSION_SEEDS[jobs$replicate[j]] + 1000L * jobs$survey[j]
  test_idx <- matched_omission(target, seed)
  fit <- held_out_fit(test_idx)
  fit$summary %>% mutate(matched_to = survey, target_n = target, replicate = jobs$replicate[j], seed = seed,
                         n_test = fit$n_test, surveys_in_test = fit$surveys_in_test, .before = 1)
}
matched <- parallel::mclapply(seq_len(nrow(jobs)), run_job, mc.cores = cores, mc.preschedule = FALSE)
failed <- vapply(matched, inherits, logical(1), what = "try-error")
if (any(failed)) stop("Matched-omission jobs failed: ", paste(which(failed), collapse = ", "))
matched <- bind_rows(matched)

# Compare each LOSO value with the matched distribution of the same size ---------
q <- function(x, p) unname(quantile(x, p, na.rm = TRUE))
matched_summary <- matched %>%
  group_by(matched_to, target_n, method, k) %>%
  summarise(
    replicates = n(),
    omitted_n_median = median(n_test),
    ari_median = median(ari, na.rm = TRUE), ari_lower = q(ari, 0.025), ari_upper = q(ari, 0.975),
    ami_median = median(ami, na.rm = TRUE), ami_lower = q(ami, 0.025), ami_upper = q(ami, 0.975),
    reference_clusters_median = median(reference_clusters),
    reference_largest_share_median = median(reference_largest_share),
    consensus_ari_median = median(som_consensus_ari, na.rm = TRUE),
    consensus_ari_lower = q(som_consensus_ari, 0.025), consensus_ari_upper = q(som_consensus_ari, 0.975),
    ratio_median = median(distance_ratio, na.rm = TRUE), ratio_lower = q(distance_ratio, 0.025), ratio_upper = q(distance_ratio, 0.975),
    .groups = "drop"
  )
comparison <- loso_summary %>%
  rename(loso_ari = ari, loso_ami = ami, loso_ratio = distance_ratio, loso_consensus_ari = som_consensus_ari) %>%
  left_join(matched_summary, by = c("held_out_campaign" = "matched_to", "method", "k")) %>%
  rowwise() %>%
  mutate(
    loso_ari_percentile = mean(matched$ari[matched$matched_to == held_out_campaign & matched$method == method & matched$k == k] <= loso_ari, na.rm = TRUE),
    loso_ami_percentile = mean(matched$ami[matched$matched_to == held_out_campaign & matched$method == method & matched$k == k] <= loso_ami, na.rm = TRUE),
    loso_ratio_percentile = mean(matched$distance_ratio[matched$matched_to == held_out_campaign & matched$method == method & matched$k == k] <= loso_ratio, na.rm = TRUE)
  ) %>%
  ungroup()

calibration <- list(loso_fits = loso_fits, loso_summary = loso_summary, matched = matched,
                    matched_summary = matched_summary, comparison = comparison)
save_rds_safe(calibration, file.path(MODEL_DIR, "r1_transfer_calibration.rds"))
write_csv_safe(comparison, file.path(TABLE_DIR, "TableS_r1_transfer_vs_matched_omission.csv"))
write_csv_safe(matched, file.path(TABLE_DIR, "TableS_r1_matched_omission_replicates.csv"))

cat("Transfer calibration completed:", nrow(jobs), "matched-omission replicates\n")
print(comparison %>% filter(k == 3L) %>%
  select(held_out_campaign, method, loso_ari, ari_median, ari_lower, ari_upper, loso_ari_percentile, loso_consensus_ari, consensus_ari_median, loso_ratio, ratio_median, ratio_upper))
