args <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "00_config.R"))

assert_packages("jsonlite")
prepared <- readRDS(file.path(RESULT_DIR, "prepared_data.rds"))
ensemble <- readRDS(file.path(MODEL_DIR, "som_ensemble_250.rds"))
resampling <- readRDS(file.path(MODEL_DIR, "cross_model_resampling_200.rds"))
transfer <- readRDS(file.path(MODEL_DIR, "campaign_transfer.rds"))
sensitivity <- readRDS(file.path(MODEL_DIR, "sensitivity_ensembles.rds"))
gradient <- readRDS(file.path(MODEL_DIR, "continuous_gradient_robustness.rds"))

methods <- c("SOM", "K-means", "Ward.D2", "GMM")
method_order <- setNames(seq_along(methods), methods)

stability_summary <- resampling$replicate_stability %>%
  group_by(method, k) %>%
  summarise(
    median_ari = median(ari), lower_ari = quantile(ari, 0.025), upper_ari = quantile(ari, 0.975),
    median_ami = median(ami), lower_ami = quantile(ami, 0.025), upper_ami = quantile(ami, 0.975),
    median_mean_jaccard = median(mean_cluster_jaccard),
    lower_mean_jaccard = quantile(mean_cluster_jaccard, 0.025),
    upper_mean_jaccard = quantile(mean_cluster_jaccard, 0.975),
    median_minimum_jaccard = median(minimum_cluster_jaccard),
    .groups = "drop"
  )

cluster_counts <- resampling$cluster_stability %>%
  dplyr::count(method, k, stability) %>%
  tidyr::complete(method, k, stability = c("stable", "uncertain", "fragile"), fill = list(n = 0L)) %>%
  tidyr::pivot_wider(names_from = stability, values_from = n, names_prefix = "clusters_")

all_stability <- stability_summary %>%
  left_join(cluster_counts, by = c("method", "k")) %>%
  left_join(resampling$internal_validity, by = c("method", "k")) %>%
  arrange(method_order[method], k)
write_csv_safe(all_stability, file.path(TABLE_DIR, "TableS_all_k_stability_and_validity.csv"))

table3 <- all_stability %>%
  filter(k %in% c(2L, 3L)) %>%
  transmute(
    Method = method, k,
    `ARI, median (95% interval)` = sprintf("%.2f (%.2f–%.2f)", median_ari, lower_ari, upper_ari),
    `AMI, median (95% interval)` = sprintf("%.2f (%.2f–%.2f)", median_ami, lower_ami, upper_ami),
    `Mean clusterwise Jaccard, median (95% interval)` = sprintf("%.2f (%.2f–%.2f)", median_mean_jaccard, lower_mean_jaccard, upper_mean_jaccard),
    `Stable / uncertain / fragile clusters` = sprintf("%d / %d / %d", clusters_stable, clusters_uncertain, clusters_fragile),
    Silhouette = round(silhouette, 2)
  )
write_csv_safe(table3, file.path(TABLE_DIR, "Table3_resampling_stability_k2_k3.csv"))

table2 <- data.frame(
  Method = methods,
  Representation = c("Hexagonal SOM codebook", "Euclidean centroids", "Ward.D2 hierarchy", "Gaussian mixture"),
  Fitting = c("Five grids × 50 seeds; 1000 online iterations", "50 starts; 200 maximum iterations", "Agglomerative Ward.D2 linkage", "Covariance model selected by BIC within training data"),
  Prediction = c("Nearest codebook unit, then codebook partition", "Nearest training centroid", "Nearest centroid of training-tree groups", "Maximum posterior probability"),
  `Shared evaluation` = rep("Same transformed variables, k=2–8 and grouped resampling indices", 4),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write_csv_safe(table2, file.path(TABLE_DIR, "Table2_controlled_model_comparison.csv"))

table4a <- transfer$summary %>%
  filter(k == 3L) %>%
  mutate(method_rank = method_order[method]) %>%
  arrange(held_out_campaign, method_rank) %>%
  transmute(
    `Held-out campaign` = held_out_campaign,
    `Test n` = n_test,
    Method = method,
    Repetitions = repetitions,
    `ARI to full-data consensus, median (95% interval)` = sprintf("%.2f (%.2f–%.2f)", median_ari, lower_ari, upper_ari),
    `AMI to full-data consensus, median` = round(median_ami, 2),
    `Test / training distance ratio, median` = round(median_distance_ratio, 2)
  )
write_csv_safe(table4a, file.path(TABLE_DIR, "Table4A_campaign_transfer_k3.csv"))

sensitivity_data <- list(
  field4 = prepared$cohorts$field4,
  nutrient4 = prepared$cohorts$nutrient4,
  nutrient_chla5 = prepared$cohorts$nutrient_chla5,
  raw_core9 = prepared$cohorts$core9,
  partial5 = prepared$cohorts$partial5
)
primary_ids <- prepared$cohorts$core9$record_id
sensitivity_agreement_rows <- list()
counter <- 0L
for (analysis_name in names(sensitivity$modules)) {
  variant_data <- sensitivity_data[[analysis_name]]
  index <- match(primary_ids, variant_data$record_id)
  if (anyNA(index)) stop("Sensitivity cohort does not contain every primary record: ", analysis_name)
  for (k_value in K_VALUES) {
    key <- paste0("k", k_value)
    counter <- counter + 1L
    sensitivity_agreement_rows[[counter]] <- data.frame(
      analysis = analysis_name,
      k = k_value,
      common_n = length(primary_ids),
      ari_to_primary_som_consensus = mclust::adjustedRandIndex(
        resampling$consensus$SOM[[key]]$labels,
        sensitivity$modules[[analysis_name]]$consensus[[key]]$labels[index]
      ),
      ami_to_primary_som_consensus = adjusted_mutual_information(
        resampling$consensus$SOM[[key]]$labels,
        sensitivity$modules[[analysis_name]]$consensus[[key]]$labels[index]
      )
    )
  }
}
sensitivity_agreement <- bind_rows(sensitivity_agreement_rows) %>%
  mutate(analysis = recode(analysis, partial5 = "partial5_native_na"))
write_csv_safe(sensitivity_agreement, file.path(TABLE_DIR, "TableS_sensitivity_to_primary_consensus.csv"))

table4b <- sensitivity$module_summary %>%
  left_join(sensitivity_agreement, by = c("analysis", "k")) %>%
  filter(k == 3L) %>%
  transmute(
    Analysis = recode(
      analysis,
      field4 = "Field chemistry (4 variables)",
      nutrient4 = "Nutrients (4 variables)",
      nutrient_chla5 = "Chl-a + nutrients (5 variables)",
      raw_core9 = "Untransformed primary matrix",
      partial5_native_na = "Native missing-distance SOM"
    ),
    n,
    `SOM models` = models,
    `ARI to ensemble consensus, median (95% interval)` = sprintf("%.2f (%.2f–%.2f)", median_ari_to_consensus, lower_ari, upper_ari),
    `ARI to primary SOM consensus` = round(ari_to_primary_som_consensus, 2),
    `Sample support, median` = round(median_sample_support, 2),
    `Assignment entropy, median` = round(median_sample_entropy, 2)
  )
write_csv_safe(table4b, file.path(TABLE_DIR, "Table4B_sensitivity_k3.csv"))

table4c <- sensitivity$exclusion %>%
  filter(k == 3L) %>%
  mutate(method_rank = method_order[method]) %>%
  arrange(method_rank) %>%
  transmute(
    Method = method, n,
    `ARI to primary consensus` = round(ari_to_primary_consensus, 2),
    `AMI to primary consensus` = round(ami_to_primary_consensus, 2)
  )
write_csv_safe(table4c, file.path(TABLE_DIR, "Table4C_source_uncertainty_exclusion_k3.csv"))

cross_summary <- resampling$cross_model_agreement %>%
  group_by(k) %>%
  summarise(
    median_ari = median(ari), minimum_ari = min(ari), maximum_ari = max(ari),
    median_ami = median(ami), minimum_ami = min(ami), maximum_ami = max(ami),
    .groups = "drop"
  )
write_csv_safe(cross_summary, file.path(TABLE_DIR, "TableS_cross_model_summary_by_k.csv"))

facts <- list(
  title = "Ensemble and cross-model validation of self-organizing maps reveals persistent gradients and fragile clusters in environmental monitoring data",
  counts = as.list(EXPECTED_COUNTS),
  campaigns_primary = as.list(table(prepared$cohorts$core9$campaign)),
  models = list(som_ensemble = 250L, grouped_resamples = 200L, som_sensitivity = 250L, campaign_transfer_rows = nrow(transfer$raw)),
  failures = list(resampling = nrow(resampling$failures), transfer = nrow(transfer$failures), sensitivity = nrow(sensitivity$failures)),
  representative_som = ensemble$representative_name,
  gradient_summary = gradient$summary,
  k3_stability = filter(all_stability, k == 3L),
  k3_cross_model = filter(resampling$cross_model_agreement, k == 3L),
  k3_transfer = filter(transfer$summary, k == 3L),
  sensitivity_to_primary = sensitivity_agreement,
  interpretation = list(
    stable_threshold = 0.75,
    uncertain_threshold = "0.60–0.75",
    fragile_threshold = 0.60,
    ecological_cluster_names_permitted = FALSE,
    reason = "No discrete partition passed joint support across grids, grouped resampling, at least two baselines, and all held-out surveys."
  )
)
jsonlite::write_json(facts, file.path(RESULT_DIR, "results_facts.json"), pretty = TRUE, auto_unbox = TRUE, digits = 8)

capture.output(sessionInfo(), file = file.path(RESULT_DIR, "sessionInfo.txt"))

cat("Main and supplementary summary tables generated\n")
cat("Ecological hard-cluster naming permitted: FALSE\n")
