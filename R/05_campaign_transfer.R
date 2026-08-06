args <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "00_config.R"))

assert_packages(c("kohonen", "mclust"))
prepared <- readRDS(file.path(RESULT_DIR, "prepared_data.rds"))
data <- prepared$cohorts$core9
resampling_file <- Sys.getenv(
  "EI_RESAMPLING_FILE",
  unset = file.path(MODEL_DIR, "cross_model_resampling_200.rds")
)
if (!file.exists(resampling_file)) stop("Resampling result is unavailable: ", resampling_file)
resampling <- readRDS(resampling_file)

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

rows <- list(); failures <- list(); counter <- 0L; failure_counter <- 0L
for (campaign_out in CAMPAIGNS) {
  train_idx <- data$campaign != campaign_out
  test_idx <- !train_idx
  train <- data[train_idx, , drop = FALSE]
  test <- data[test_idx, , drop = FALSE]
  transform <- transform_fit(train)
  x_train <- transform_apply(train, transform)
  x_test <- transform_apply(test, transform)
  x_combined <- rbind(x_train, x_test)

  for (repeat_id in seq_along(TRANSFER_SEEDS)) {
    seed <- TRANSFER_SEEDS[repeat_id]
    g <- 1L + ((repeat_id - 1L) %% nrow(GRIDS))
    set.seed(seed)
    som <- tryCatch(
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
    for (k in K_VALUES) {
      key <- paste0("k", k)
      reference <- resampling$consensus$SOM[[key]]$labels[test_idx]
      if (inherits(som, "error")) {
        failure_counter <- failure_counter + 1L
        failures[[failure_counter]] <- data.frame(held_out_campaign = campaign_out, method = "SOM", repeat_id = repeat_id, k = k, message = conditionMessage(som))
      } else {
        labels <- partition_som(som, k, x_combined)
        train_labels <- labels[seq_len(nrow(x_train))]
        test_labels <- labels[nrow(x_train) + seq_len(nrow(x_test))]
        counter <- counter + 1L
        rows[[counter]] <- data.frame(
          held_out_campaign = campaign_out, n_train = nrow(train), n_test = nrow(test),
          method = "SOM", repeat_id = repeat_id, seed = seed,
          configuration = sprintf("%dx%d", GRIDS$xdim[g], GRIDS$ydim[g]), k = k,
          ari_to_full_consensus = mclust::adjustedRandIndex(reference, test_labels),
          ami_to_full_consensus = adjusted_mutual_information(reference, test_labels),
          distance_ratio = centroid_distance_ratio(x_train, x_test, train_labels, test_labels)
        )
      }

      km <- fit_kmeans_predict(x_train, x_combined, k, seed + k)
      km_train <- km$labels[seq_len(nrow(x_train))]
      km_test <- km$labels[nrow(x_train) + seq_len(nrow(x_test))]
      km_reference <- resampling$consensus[["K-means"]][[key]]$labels[test_idx]
      counter <- counter + 1L
      rows[[counter]] <- data.frame(
        held_out_campaign = campaign_out, n_train = nrow(train), n_test = nrow(test),
        method = "K-means", repeat_id = repeat_id, seed = seed,
        configuration = "50 starts", k = k,
        ari_to_full_consensus = mclust::adjustedRandIndex(km_reference, km_test),
        ami_to_full_consensus = adjusted_mutual_information(km_reference, km_test),
        distance_ratio = centroid_distance_ratio(x_train, x_test, km_train, km_test)
      )
    }
  }

  for (k in K_VALUES) {
    key <- paste0("k", k)
    ward <- fit_ward_predict(x_train, x_combined, k)
    ward_train <- ward$labels[seq_len(nrow(x_train))]
    ward_test <- ward$labels[nrow(x_train) + seq_len(nrow(x_test))]
    ward_reference <- resampling$consensus[["Ward.D2"]][[key]]$labels[test_idx]
    counter <- counter + 1L
    rows[[counter]] <- data.frame(
      held_out_campaign = campaign_out, n_train = nrow(train), n_test = nrow(test),
      method = "Ward.D2", repeat_id = 1L, seed = NA_integer_, configuration = "Ward linkage", k = k,
      ari_to_full_consensus = mclust::adjustedRandIndex(ward_reference, ward_test),
      ami_to_full_consensus = adjusted_mutual_information(ward_reference, ward_test),
      distance_ratio = centroid_distance_ratio(x_train, x_test, ward_train, ward_test)
    )

    gmm <- fit_gmm_predict(x_train, x_combined, k)
    if (!is.na(gmm$failure)) {
      failure_counter <- failure_counter + 1L
      failures[[failure_counter]] <- data.frame(held_out_campaign = campaign_out, method = "GMM", repeat_id = 1L, k = k, message = gmm$failure)
    } else {
      gmm_train <- gmm$labels[seq_len(nrow(x_train))]
      gmm_test <- gmm$labels[nrow(x_train) + seq_len(nrow(x_test))]
      gmm_reference <- resampling$consensus$GMM[[key]]$labels[test_idx]
      counter <- counter + 1L
      rows[[counter]] <- data.frame(
        held_out_campaign = campaign_out, n_train = nrow(train), n_test = nrow(test),
        method = "GMM", repeat_id = 1L, seed = NA_integer_, configuration = gmm$model_name, k = k,
        ari_to_full_consensus = mclust::adjustedRandIndex(gmm_reference, gmm_test),
        ami_to_full_consensus = adjusted_mutual_information(gmm_reference, gmm_test),
        distance_ratio = centroid_distance_ratio(x_train, x_test, gmm_train, gmm_test)
      )
    }
  }
  message("Campaign transfer completed: ", campaign_out)
}

transfer <- bind_rows(rows)
failure_table <- if (length(failures)) bind_rows(failures) else data.frame(held_out_campaign = character(), method = character(), repeat_id = integer(), k = integer(), message = character())
summary <- transfer %>%
  group_by(held_out_campaign, n_test, method, k) %>%
  summarise(
    repetitions = n(),
    median_ari = median(ari_to_full_consensus),
    lower_ari = quantile(ari_to_full_consensus, 0.025),
    upper_ari = quantile(ari_to_full_consensus, 0.975),
    median_ami = median(ami_to_full_consensus),
    median_distance_ratio = median(distance_ratio),
    .groups = "drop"
  )

save_rds_safe(list(raw = transfer, summary = summary, failures = failure_table), file.path(MODEL_DIR, "campaign_transfer.rds"))
write_csv_safe(transfer, file.path(TABLE_DIR, "TableS_campaign_transfer_all.csv"))
write_csv_safe(summary, file.path(TABLE_DIR, "Table4_campaign_transfer_summary.csv"))
write_csv_safe(failure_table, file.path(TABLE_DIR, "TableS_campaign_transfer_failures.csv"))

cat("Campaign transfer rows:", nrow(transfer), "\n")
cat("Campaign transfer failures:", nrow(failure_table), "\n")
