args <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "00_config.R"))

# Display tables and a numbers file for the first-round revision. Supplementary
# Tables S15-S21 are new; Tables S3 and S11 and main Table 4 are revised. Every
# number quoted in the revised text is collected in results_facts_r1.json so
# that the manuscript can be checked against the analysis automatically.

assert_packages(c("jsonlite", "kohonen"))
prepared <- readRDS(file.path(RESULT_DIR, "prepared_data.rds"))
ensemble <- readRDS(file.path(MODEL_DIR, "som_ensemble_250.rds"))
resampling <- readRDS(file.path(MODEL_DIR, "cross_model_resampling_200.rds"))
gradient <- readRDS(file.path(MODEL_DIR, "continuous_gradient_robustness.rds"))
topology <- readRDS(file.path(MODEL_DIR, "r1_som_topology.rds"))
stages <- readRDS(file.path(MODEL_DIR, "r1_partition_stages.rds"))
metrics <- readRDS(file.path(MODEL_DIR, "r1_metric_checks.rds"))
calibration <- readRDS(file.path(MODEL_DIR, "r1_transfer_calibration.rds"))
like <- readRDS(file.path(MODEL_DIR, "r1_sensitivity_like_for_like.rds"))
methods <- c("SOM", "K-means", "Ward.D2", "GMM")
method_rank <- function(x) match(x, methods)

f3 <- function(x) formatC(x, format = "f", digits = 3)
f2 <- function(x) formatC(x, format = "f", digits = 2)
interval3 <- function(m, lo, hi) sprintf("%s (%s–%s)", f3(m), f3(lo), f3(hi))
pct <- function(x) sprintf("%.0f%%", 100 * x)

# Table S3 (revised): quantisation error as reported and as Euclidean distance ----
x <- prepared$x
euclidean_qe <- vapply(ensemble$models, function(model) {
  mean(nearest_units(x, model$codes[[1]])$distance)
}, numeric(1))
table_s3 <- ensemble$metrics %>%
  mutate(euclidean_qe = euclidean_qe, grid = factor(grid, levels = c("5x5", "7x5", "8x6", "10x7", "12x8"))) %>%
  group_by(grid) %>%
  summarise(
    Models = n(),
    `Quantisation error, mean squared distance` = f3(median(quantization_error)),
    `Quantisation error, mean Euclidean distance` = f3(median(euclidean_qe)),
    `Topographic error` = f3(median(topographic_error)),
    `10-neighbour preservation` = f3(median(neighbour_preservation_10)),
    `Empty nodes` = median(empty_nodes),
    .groups = "drop"
  ) %>%
  rename(`SOM grid` = grid)
write_csv_safe(table_s3, file.path(TABLE_DIR, "TableS3_som_grid_summary_r1_display.csv"))

# Table S15: sample-level support, margin and entropy ---------------------------
es <- metrics$entropy_summary
max_self_shift <- max(abs(c(es$median_entropy - es$median_entropy_excluding_self,
                            es$median_support - es$median_support_excluding_self,
                            es$median_margin - es$median_margin_excluding_self)))
table_s15 <- es %>%
  arrange(method_rank(method), k) %>%
  transmute(
    Method = method, k,
    `Normalised entropy, median (IQR)` = sprintf("%s (%s–%s)", f3(median_entropy), f3(q25_entropy), f3(q75_entropy)),
    `Normalised entropy, mean` = f3(mean_entropy),
    `Support, median` = f3(median_support),
    `Margin, median` = f3(median_margin),
    `Singleton consensus clusters` = singleton_events
  )
write_csv_safe(table_s15, file.path(TABLE_DIR, "TableS15_sample_metrics_display.csv"))

# Table S16: clusterwise stability classes from median and mean Jaccard -----------
table_s16 <- metrics$jaccard_classes %>%
  arrange(method_rank(method), k) %>%
  transmute(
    Method = method, k,
    `Stable / uncertain / fragile, median Jaccard` = classes_median,
    `Stable / uncertain / fragile, mean Jaccard` = classes_mean,
    `Clusters reclassified` = clusters_changed,
    `Stable / uncertain / fragile, median Jaccard on training events only` = classes_training_only,
    `Clusters reclassified (training events only)` = clusters_changed_training_only
  )
write_csv_safe(table_s16, file.path(TABLE_DIR, "TableS16_jaccard_median_vs_mean_display.csv"))

# Table S17: reproducibility of the SOM map layout ------------------------------
table_s17 <- topology$summary %>%
  mutate(order = match(design, c("Full-data SOM ensemble", "Grouped-resample SOMs", "Grouped-resample PCA (two axes)",
                                 "Grouped-resample PCA quantised to the SOM grid"))) %>%
  arrange(order, grid != "All grids", grid, comparison, metric) %>%
  transmute(
    Fits = recode(design,
      "Full-data SOM ensemble" = "250 full-data SOMs",
      "Grouped-resample SOMs" = "200 grouped-resample SOMs",
      "Grouped-resample PCA (two axes)" = "200 grouped-resample PCAs (two axes)",
      "Grouped-resample PCA quantised to the SOM grid" = "200 grouped-resample PCAs quantised to the matching SOM grid"),
    Pairs = paste0(comparison, " (", pairs, ")"),
    Grid = grid,
    Metric = recode(metric, "Map-distance concordance" = "Map-layout concordance"),
    `Median (95% interval)` = interval3(median, lower_95, upper_95)
  )
write_csv_safe(table_s17, file.path(TABLE_DIR, "TableS17_som_map_reproducibility_display.csv"))

# Table S18: SOM second-level clustering and consensus level --------------------
table_s18 <- stages$som_second_level %>%
  mutate(design_order = ifelse(grepl("^Grouped", design), 1L, 2L),
         variant_order = ifelse(grepl("^Unweighted", variant), 1L, 2L)) %>%
  arrange(design_order, k, variant_order) %>%
  transmute(
    Design = ifelse(grepl("^Grouped", design), "200 grouped resamples", "250 grid–seed maps"),
    `Prototype clustering` = ifelse(grepl("^Unweighted", variant), "Unweighted, all prototypes (primary)", "Hit-weighted, occupied prototypes"),
    k,
    `ARI to own consensus, median (95% interval)` = interval3(median_ari, lower_ari, upper_ari),
    `Stable / uncertain / fragile` = classes_median_jaccard,
    `Consensus cluster sizes` = consensus_sizes,
    `Fits with a cluster of ≤5 events` = pct(share_fits_smallest_le5),
    `ARI with primary consensus` = f3(ari_to_primary_consensus)
  )
write_csv_safe(table_s18, file.path(TABLE_DIR, "TableS18_som_second_level_display.csv"))
write_csv_safe(stages$consensus_level %>% mutate(across(-k, f3)), file.path(TABLE_DIR, "TableS18b_consensus_level_display.csv"))

# Table S19: Ward.D2 and K-means allocation of training events -------------------
table_s19 <- stages$ward_allocation %>%
  transmute(
    k,
    `Tree vs nearest-centroid labels, median ARI (minimum)` = sprintf("%s (%s)", f3(ward_tree_vs_nearest_median_ari), f3(ward_tree_vs_nearest_min_ari)),
    `Identical memberships, median (minimum)` = sprintf("%s (%s)", pct(ward_tree_vs_nearest_median_identical), pct(ward_tree_vs_nearest_min_identical)),
    `Ward.D2 with tree memberships: ARI to own consensus, median (95% interval)` = interval3(median_ari, lower_ari, upper_ari),
    `Stable / uncertain / fragile` = classes_median_jaccard,
    `ARI with primary Ward.D2 consensus` = f3(ari_tree_consensus_vs_primary_ward),
    `ARI with K-means consensus` = f3(ari_tree_consensus_vs_kmeans),
    `K-means fits with identical memberships` = pct(kmeans_share_fits_identical)
  )
write_csv_safe(table_s19, file.path(TABLE_DIR, "TableS19_ward_kmeans_allocation_display.csv"))

# Table S20: survey transfer against matched-size random omission ----------------
# A held-out reference with a single cluster gives an ARI of 0 or 1 only, so such
# combinations are flagged per method and k and left out of the counts.
comparison <- calibration$comparison %>% mutate(informative = reference_clusters > 1L)
below_counts <- comparison %>%
  group_by(held_out_campaign, method) %>%
  summarise(k_informative = sum(informative),
            k_below = sum(informative & loso_ari < ari_lower, na.rm = TRUE),
            k_below_ami = sum(informative & loso_ami < ami_lower, na.rm = TRUE),
            k_ratio_above = sum(loso_ratio > ratio_upper, na.rm = TRUE), .groups = "drop")
table_s20 <- comparison %>%
  filter(k == 3L) %>%
  left_join(below_counts, by = c("held_out_campaign", "method")) %>%
  arrange(held_out_campaign, method_rank(method)) %>%
  transmute(
    `Held-out survey` = held_out_campaign,
    `Held-out n` = n_test,
    Method = method,
    `Reference clusters (largest share): LOSO; matched median` = sprintf("%d (%s); %s (%s)", reference_clusters,
      pct(reference_largest_share), formatC(reference_clusters_median, format = "fg"), pct(reference_largest_share_median)),
    `ARI: LOSO (percentile); matched median (95% interval)` = sprintf("%s; %s",
      ifelse(informative, sprintf("%s (%s)", f3(loso_ari), pct(loso_ari_percentile)), sprintf("%s (single-cluster reference)", f3(loso_ari))),
      interval3(ari_median, ari_lower, ari_upper)),
    `AMI: LOSO (percentile); matched median (95% interval)` = sprintf("%s; %s",
      ifelse(informative, sprintf("%s (%s)", f3(loso_ami), pct(loso_ami_percentile)), sprintf("%s (single-cluster reference)", f3(loso_ami))),
      interval3(ami_median, ami_lower, ami_upper)),
    `SOM 10-map consensus ARI: LOSO; matched median (95% interval)` = ifelse(method == "SOM",
      sprintf("%s; %s", f3(loso_consensus_ari), interval3(consensus_ari_median, consensus_ari_lower, consensus_ari_upper)), ""),
    `Distance ratio: LOSO; matched median (95% interval)` = sprintf("%s; %s", f2(loso_ratio), sprintf("%s (%s–%s)", f2(ratio_median), f2(ratio_lower), f2(ratio_upper))),
    `k values below matched range, ARI / AMI (of k with a multi-cluster reference)` = ifelse(k_informative == 0L, "–",
      sprintf("%d / %d (of %d)", k_below, k_below_ami, k_informative))
  )
write_csv_safe(table_s20, file.path(TABLE_DIR, "TableS20_transfer_matched_omission_display.csv"))

# Table S21: variable modules with variables and events changed one at a time ----
table_s21 <- like$modules %>%
  filter(k == 3L) %>%
  transmute(
    Analysis = analysis, Events = cohort, n,
    `Within-variant ARI, median (95% interval)` = interval3(within_median, within_lower, within_upper),
    `ARI with primary grid–seed consensus` = f3(agreement_primary_grid_seed),
    `ARI with primary resampling consensus` = f3(agreement_primary_resampling),
    `ARI between common-event and full-cohort module consensus` = ifelse(is.na(agreement_with_full_module_cohort) | agreement_with_full_module_cohort == 1, "", f3(agreement_with_full_module_cohort))
  )
write_csv_safe(table_s21, file.path(TABLE_DIR, "TableS21_modules_like_for_like_display.csv"))
write_csv_safe(like$modules, file.path(TABLE_DIR, "TableS21_modules_like_for_like_all_k.csv"))

# Table S11 (revised): exclusion of the 2024-2025 survey -------------------------
table_s11 <- like$exclusion %>%
  arrange(method_rank(method), k) %>%
  transmute(
    Method = method, k,
    `ARI with primary consensus (374 events)` = f3(ensemble_consensus_ari),
    `Random omission of ~16 events, median (95% interval)` = interval3(random_omission_median, random_omission_lower, random_omission_upper),
    `SOM single maps, median (95% interval)` = ifelse(method == "SOM", interval3(single_map_median, single_map_lower, single_map_upper), ""),
    `Single map reported in submitted version` = ifelse(method == "SOM", f3(submitted_single_map_ari), "")
  )
write_csv_safe(table_s11, file.path(TABLE_DIR, "TableS11_survey_exclusion_r1_display.csv"))

# Main Table 4 (revised): diagnostic heuristic with explicit evidence guides -----
# The rows follow the joint evidence rule of Section 2.2 and are applied in order
# to the SOM consensus at each k, so that every pattern falls in exactly one row.
table4 <- data.frame(
  `Evidence pattern at a given k (interpretive guide used here)` = c(
    "Dominant geometry is not reproducible (median PCA score-distance concordance between grouped resamples < 0.95)",
    "Geometry is reproducible, but any SOM cluster is fragile (median Jaccard < 0.60) or two or more are uncertain (0.60–0.75)",
    "No fragile and at most one uncertain SOM cluster, but the joint evidence rule is not met (one uncertain cluster, SOM consensus ARI < 0.80 with a contrasting method, or held-out agreement below the matched 95% range for an informative survey)",
    "Joint evidence rule met: every SOM cluster ≥ 0.75, SOM consensus ARI ≥ 0.80 with every contrasting method, and held-out agreement not below the matched 95% range for any informative survey"
  ),
  Interpretation = c(
    "Representation is unstable",
    "Persistent gradients with unsupported boundaries",
    "Conditional segmentation",
    "Candidate ecological classes"
  ),
  `Defensible output` = c(
    "No environmental classification",
    "Continuous scores, best matching units on a stated reference map, support and margins",
    "Numerical groups with sensitivity bounds",
    "Provisionally profiled groups within the fitted archive"
  ),
  `Monitoring or analysis action` = c(
    "Review variables, transformations, sampling blocks and SOM settings before interpretation",
    "Treat small margins as possible transitions; resample only where they coincide with sparse coverage",
    "Limit inference to the fitted domain and disclose the cluster concept",
    "Validate with new surveys or external variables before operational naming"
  ),
  check.names = FALSE
)
write_csv_safe(table4, file.path(TABLE_DIR, "Table4_decision_heuristic_r1.csv"))

# Numbers quoted in the revised manuscript --------------------------------------
pick <- function(df, ...) { out <- dplyr::filter(df, ...); if (nrow(out) != 1) stop("Expected one row"); out }
tsum <- topology$summary
som2 <- stages$som_second_level
ward <- stages$ward_allocation
cl <- stages$consensus_level
exc <- like$exclusion
mods <- like$modules
k3 <- comparison %>% filter(k == 3L)
facts <- list(
  euclidean_qe_5x5 = median(euclidean_qe[ensemble$metrics$grid == "5x5"]),
  euclidean_qe_12x8 = median(euclidean_qe[ensemble$metrics$grid == "12x8"]),
  map_concordance_full = pick(tsum, design == "Full-data SOM ensemble", comparison == "All pairs", grid == "All grids")[c("median", "lower_95", "upper_95")],
  map_concordance_full_same_grid = pick(tsum, design == "Full-data SOM ensemble", comparison == "Same-grid pairs", grid == "All grids", metric == "Map-distance concordance")[c("median", "lower_95", "upper_95")],
  map_concordance_resample = pick(tsum, design == "Grouped-resample SOMs", comparison == "All pairs", grid == "All grids")[c("median", "lower_95", "upper_95")],
  neighbour_overlap_full = pick(tsum, design == "Full-data SOM ensemble", grid == "All grids", metric == "Map-neighbour overlap")[c("median", "lower_95", "upper_95")],
  neighbour_overlap_resample = pick(tsum, design == "Grouped-resample SOMs", grid == "All grids", metric == "Map-neighbour overlap")[c("median", "lower_95", "upper_95")],
  neighbour_overlap_chance = pick(tsum, design == "Full-data SOM ensemble", grid == "All grids", metric == "Map-neighbour overlap, chance")[c("median")],
  pca_between_resample = pick(tsum, design == "Grouped-resample PCA (two axes)")[c("median", "lower_95", "upper_95")],
  som_k3_resample_weighted = pick(som2, grepl("^Grouped", design), grepl("^Hit", variant), k == 3L)[c("median_ari", "lower_ari", "upper_ari", "classes_median_jaccard", "ari_to_primary_consensus")],
  som_k3_grid_seed_unweighted = pick(som2, grepl("^Grid", design), grepl("^Unweighted", variant), k == 3L)[c("median_ari", "lower_ari", "upper_ari")],
  som_k3_grid_seed_weighted = pick(som2, grepl("^Grid", design), grepl("^Hit", variant), k == 3L)[c("median_ari", "lower_ari", "upper_ari")],
  som_small_cluster_share = som2 %>% filter(k >= 4L) %>% group_by(design, variant) %>%
    summarise(min = min(share_fits_smallest_le5), max = max(share_fits_smallest_le5), .groups = "drop"),
  som_weighted_fragile_by_k = som2 %>% filter(grepl("^Grouped", design), grepl("^Hit", variant)) %>% select(k, fragile_clusters, consensus_sizes),
  consensus_level = cl,
  ward_k3 = pick(ward, k == 3L)[c("median_ari", "lower_ari", "upper_ari", "classes_median_jaccard", "ward_tree_vs_nearest_median_ari", "ward_tree_vs_nearest_median_identical", "ward_tree_vs_nearest_min_identical", "ari_tree_consensus_vs_primary_ward", "ari_tree_consensus_vs_kmeans")],
  ward_identical_range = range(ward$ward_tree_vs_nearest_median_identical),
  ward_min_identical_range = range(ward$ward_tree_vs_nearest_min_identical),
  jaccard_classes_changed = metrics$jaccard_classes %>% filter(clusters_changed > 0) %>% select(method, k, classes_median, classes_mean),
  max_self_pair_shift = max_self_shift,
  entropy_k3 = metrics$entropy_summary %>% filter(k == 3L) %>% select(method, median_entropy, q25_entropy, q75_entropy, mean_entropy),
  entropy_order = metrics$entropy_rank,
  transfer_k3 = k3 %>% select(held_out_campaign, method, loso_ari, ari_median, ari_lower, ari_upper, loso_ari_percentile,
                              loso_ami, ami_median, ami_lower, ami_upper, loso_ami_percentile,
                              reference_clusters, reference_largest_share, reference_clusters_median, reference_largest_share_median,
                              loso_consensus_ari, consensus_ari_median, consensus_ari_lower, consensus_ari_upper,
                              loso_ratio, ratio_median, ratio_lower, ratio_upper),
  transfer_below_counts = below_counts,
  transfer_tallies = list(
    informative_surveys = comparison %>% filter(held_out_campaign != "2024-2025", informative) %>%
      summarise(combinations = n(), ari_below = sum(loso_ari < ari_lower), ami_below = sum(loso_ami < ami_lower)),
    all_surveys = comparison %>% filter(informative) %>%
      summarise(combinations = n(), ari_below = sum(loso_ari < ari_lower), ami_below = sum(loso_ami < ami_lower)),
    single_cluster_references = comparison %>% filter(!informative) %>% select(held_out_campaign, method, k)
  ),
  matched_ratio_range_k3 = range(k3$ratio_median),
  matched_ratio_range_k3_informative = range(k3$ratio_median[k3$held_out_campaign != "2024-2025"]),
  map_concordance_resample_same_grid = pick(tsum, design == "Grouped-resample SOMs", comparison == "Same-grid pairs", grid == "All grids", metric == "Map-distance concordance")[c("median", "lower_95", "upper_95")],
  quantised_pca = pick(tsum, design == "Grouped-resample PCA quantised to the SOM grid", comparison == "All pairs")[c("median", "lower_95", "upper_95")],
  quantised_pca_same_grid = pick(tsum, design == "Grouped-resample PCA quantised to the SOM grid", comparison == "Same-grid pairs")[c("median", "lower_95", "upper_95")],
  component_relation_resample = pick(tsum, metric == "Component-relation correlation")[c("median", "lower_95", "upper_95")],
  jaccard_training_only = metrics$jaccard_classes %>% select(method, k, classes_median, classes_training_only, clusters_changed_training_only, fragile_training_only),
  som_stability_by_k = som2 %>% filter(grepl("^Grouped", design)) %>% select(variant, k, median_ari, classes_median_jaccard, ari_to_primary_consensus),
  ward_tree_by_k = ward %>% select(k, median_ari, classes_median_jaccard),
  exclusion = exc %>% filter(k %in% 2:8) %>% select(method, k, ensemble_consensus_ari, submitted_single_map_ari,
                                                     single_map_median, single_map_lower, single_map_upper,
                                                     random_omission_median, random_omission_lower, random_omission_upper),
  modules_k3 = mods %>% filter(k == 3L),
  missing_identity = like$missing_identity
)
jsonlite::write_json(facts, file.path(RESULT_DIR, "results_facts_r1.json"), auto_unbox = TRUE, pretty = TRUE, digits = 8)
cat("R1 display tables and facts written\n")
