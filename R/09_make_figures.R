# Recreate manuscript Figs. 3-6 and the analytical supplementary figures.
args <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "00_config.R"))

assert_packages(c("patchwork", "svglite", "scales"))
library(patchwork)
prepared <- readRDS(file.path(RESULT_DIR, "prepared_data.rds"))
ensemble <- readRDS(file.path(MODEL_DIR, "som_ensemble_250.rds"))
resampling <- readRDS(file.path(MODEL_DIR, "cross_model_resampling_200.rds"))
transfer <- readRDS(file.path(MODEL_DIR, "campaign_transfer.rds"))
sensitivity <- readRDS(file.path(MODEL_DIR, "sensitivity_ensembles.rds"))
gradient <- readRDS(file.path(MODEL_DIR, "continuous_gradient_robustness.rds"))

# Use colour-blind-friendly palettes with line and symbol redundancy.
method_colours <- c("SOM" = "#0072B2", "K-means" = "#E69F00", "Ward.D2" = "#009E73", "GMM" = "#CC79A7")
campaign_colours <- c("2012-2013" = "#0072B2", "2017-2018" = "#E69F00", "2020-2021" = "#009E73", "2024-2025" = "#D55E00")
campaign_linetypes <- c("2012-2013" = "solid", "2017-2018" = "dashed", "2020-2021" = "dotdash", "2024-2025" = "dotted")
campaign_shapes <- c("2012-2013" = 16, "2017-2018" = 17, "2020-2021" = 15, "2024-2025" = 18)
stability_colours <- c("stable" = "#009E73", "uncertain" = "#E69F00", "fragile" = "#D55E00")
base_theme <- theme_bw(base_family = "Helvetica", base_size = 8) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#E8E8E8", linewidth = 0.25),
    strip.background = element_rect(fill = "#F2F2F2", colour = "#B5B5B5", linewidth = 0.3),
    strip.text = element_text(face = "bold"),
    legend.key.height = unit(3.5, "mm"),
    plot.title = element_text(face = "bold", size = 9),
    plot.tag = element_text(face = "bold", size = 10)
  )

# Export each analytical figure in vector, print and review formats.
save_multiformat <- function(plot, stem, width_mm = 180, height_mm = 120, directory = FIG_MAIN_DIR) {
  ggsave(file.path(directory, paste0(stem, ".pdf")), plot, width = width_mm, height = height_mm,
         units = "mm", device = grDevices::pdf, useDingbats = FALSE)
  ggsave(file.path(directory, paste0(stem, ".svg")), plot, width = width_mm, height = height_mm,
         units = "mm", device = svglite::svglite)
  ggsave(file.path(directory, paste0(stem, ".png")), plot, width = width_mm, height = height_mm,
         units = "mm", dpi = 600, bg = "white")
  ggsave(file.path(directory, paste0(stem, ".tiff")), plot, width = width_mm, height = height_mm,
         units = "mm", dpi = 600, compression = "lzw", bg = "white")
}

write_caption <- function(stem, text_value, directory = FIG_MAIN_DIR) {
  writeLines(text_value, file.path(directory, paste0(stem, "_caption.md")), useBytes = TRUE)
}

# Fig. 3: SOM component planes and topology diagnostics ----------------------
representative <- ensemble$models[[ensemble$representative_index]]
centres <- as.data.frame(representative$grid$pts)
codes <- as.data.frame(representative$codes[[1]])
var_labels <- c(
  pH = "pH", WT_C = "Water temperature", DO_mgL = "Dissolved oxygen",
  Cond_uScm = "Conductivity", Chla_ugL = "log10 Chl-a", TN_mgL = "log10 TN",
  NO3_mgL = "log10 NO3-N", TP_mgL = "log10 TP", SRP_mgL = "log10 SRP"
)
component <- bind_cols(centres, codes) %>%
  pivot_longer(cols = all_of(VARS), names_to = "variable", values_to = "code") %>%
  mutate(variable = factor(var_labels[variable], levels = unname(var_labels)))

hex_polygons <- function(component_data, radius = 0.56) {
  angles <- seq(0, 2 * pi, length.out = 7) + pi / 6
  do.call(rbind, lapply(seq_len(nrow(component_data)), function(i) {
    data.frame(
      variable = component_data$variable[i], node = i,
      x = component_data$x[i] + radius * cos(angles),
      y = component_data$y[i] + radius * sin(angles),
      code = component_data$code[i]
    )
  }))
}
component_polygons <- component %>%
  group_split(variable) %>%
  lapply(hex_polygons) %>% bind_rows()

p_components <- ggplot(component_polygons, aes(x, y, group = interaction(variable, node), fill = code)) +
  geom_polygon(colour = "white", linewidth = 0.08) +
  facet_wrap(~variable, ncol = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
                       limits = c(-3, 3), oob = scales::squish,
                       name = "Prototype\n(z score)") +
  coord_equal() +
  labs(subtitle = "Representative 12 × 8 SOM") +
  theme_void(base_family = "Helvetica", base_size = 8) +
  theme(
    strip.text = element_text(face = "bold", margin = margin(2, 1, 2, 1)),
    legend.position = "right", plot.subtitle = element_text(size = 8, face = "italic"),
    panel.spacing = unit(1.5, "mm")
  )

topology_long <- ensemble$metrics %>%
  select(grid, quantization_error, topographic_error, neighbour_preservation_10) %>%
  pivot_longer(-grid, names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
    quantization_error = "Quantization error",
    topographic_error = "Topographic error",
    neighbour_preservation_10 = "10-neighbour preservation"
  ), grid = factor(grid, levels = c("5x5", "7x5", "8x6", "10x7", "12x8")))
p_topology <- ggplot(topology_long, aes(grid, value, colour = grid)) +
  geom_violin(aes(fill = grid), alpha = 0.14, colour = NA, scale = "width") +
  geom_boxplot(width = 0.24, outlier.shape = NA, fill = "white", linewidth = 0.35) +
  geom_point(position = position_jitter(width = 0.10, height = 0), size = 0.45, alpha = 0.35) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_colour_viridis_d(option = "C", end = 0.88) +
  scale_fill_viridis_d(option = "C", end = 0.88) +
  labs(x = "SOM grid", y = NULL) + base_theme +
  theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1))

fig3 <- p_components / p_topology + plot_layout(heights = c(2.5, 1)) + plot_annotation(tag_levels = "A")
save_multiformat(fig3, "Fig3_som_components_and_topology", 180, 194)
write_caption("Fig3_som_components_and_topology", paste(
  "SOM component planes and ensemble topology diagnostics.",
  "(A) Prototype values for the representative model chosen only from quantization error, topographic error and neighbourhood preservation; colours show standardised transformed values clipped at ±3 for display, not water-quality classes.",
  "(B) Distributions across 50 seeds for each of five prespecified hexagonal grids."
))

# Fig. 4: stability over k and across model families --------------------------
stability_summary <- resampling$replicate_stability %>%
  group_by(method, k) %>%
  summarise(
    median = median(mean_cluster_jaccard),
    lower = quantile(mean_cluster_jaccard, 0.025), upper = quantile(mean_cluster_jaccard, 0.975),
    .groups = "drop"
  ) %>% mutate(method = factor(method, levels = names(method_colours)))
p4a <- ggplot(stability_summary, aes(k, median, colour = method, fill = method)) +
  geom_hline(yintercept = c(0.60, 0.75), linetype = c("dotted", "dashed"), colour = "#555555", linewidth = 0.35) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.75) + geom_point(size = 1.6) +
  scale_colour_manual(values = method_colours) + scale_fill_manual(values = method_colours) +
  scale_x_continuous(breaks = K_VALUES) + coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Number of clusters (k)", y = "Mean clusterwise Jaccard", colour = NULL, fill = NULL) + base_theme +
  theme(legend.position = "top")

cluster_counts <- resampling$cluster_stability %>%
  dplyr::count(method, k, stability) %>%
  tidyr::complete(method, k, stability = c("stable", "uncertain", "fragile"), fill = list(n = 0L)) %>%
  group_by(method, k) %>%
  mutate(total = sum(n), fraction = n / total) %>% ungroup() %>%
  mutate(method = factor(method, levels = rev(names(method_colours))))
p4b <- ggplot(cluster_counts, aes(k, method, fill = stability, alpha = fraction)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = n), size = 2.4) +
  facet_grid(stability ~ ., switch = "y") +
  scale_fill_manual(values = stability_colours) +
  scale_alpha_continuous(range = c(0.12, 1), guide = "none") +
  scale_x_continuous(breaks = K_VALUES) +
  labs(x = "Number of clusters (k)", y = NULL, fill = NULL) + base_theme +
  theme(legend.position = "none", strip.placement = "outside", strip.background = element_blank())

cross <- resampling$cross_model_agreement %>%
  mutate(pair = paste(method_1, method_2, sep = " - "))
p4c <- ggplot(cross, aes(k, ari, colour = pair, linetype = pair)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "#777777") +
  geom_line(linewidth = 0.62) + geom_point(size = 1.25) +
  scale_x_continuous(breaks = K_VALUES) +
  scale_colour_viridis_d(option = "D", end = 0.90) +
  labs(x = "Number of clusters (k)", y = "Cross-model ARI", colour = NULL, linetype = NULL) + base_theme +
  theme(legend.position = "right", legend.text = element_text(size = 6.5))

entropy_summary <- resampling$sample_metrics %>%
  group_by(method, k) %>%
  summarise(median = median(entropy), lower = quantile(entropy, 0.25), upper = quantile(entropy, 0.75), .groups = "drop") %>%
  mutate(method = factor(method, levels = names(method_colours)))
p4d <- ggplot(entropy_summary, aes(k, median, colour = method, fill = method)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.75) + geom_point(size = 1.6) +
  scale_colour_manual(values = method_colours) + scale_fill_manual(values = method_colours) +
  scale_x_continuous(breaks = K_VALUES) + coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Number of clusters (k)", y = "Normalised assignment entropy", colour = NULL, fill = NULL) + base_theme +
  theme(legend.position = "none")

fig4 <- (p4a | p4b) / (p4c | p4d) + plot_annotation(tag_levels = "A")
save_multiformat(fig4, "Fig4_resampling_and_cross_model_stability", 180, 158)
write_caption("Fig4_resampling_and_cross_model_stability", paste(
  "Hard-partition stability and cross-model agreement across k = 2–8.",
  "(A) Median and 95% empirical intervals of mean clusterwise Jaccard similarity over 200 grouped resamples; horizontal lines mark the prespecified fragile/uncertain and uncertain/stable boundaries.",
  "(B) Numbers of stable, uncertain and fragile clusters.",
  "(C) Pairwise adjusted Rand index (ARI) between full resampling consensus partitions.",
  "(D) Median and interquartile range of sample assignment entropy. Internal indices are not used as pass/fail criteria."
))

# Fig. 5: continuous geometry and sample-level uncertainty -------------------
scores <- as.data.frame(gradient$full_pca$x[, 1:2])
names(scores) <- c("PC1", "PC2")
scores <- bind_cols(prepared$cohorts$core9 %>% select(campaign), scores)
som_k3 <- resampling$sample_metrics %>% filter(method == "SOM", k == 3L) %>% arrange(row_index)
scores <- bind_cols(scores, som_k3 %>% select(support, entropy, consensus_cluster))
loadings <- as.data.frame(gradient$full_pca$rotation[, 1:2])
names(loadings) <- c("PC1", "PC2")
loadings$variable <- rownames(loadings)
loadings$label <- c("pH", "WT", "DO", "Cond", "Chl-a", "TN", "NO3-N", "TP", "SRP")
loadings$label_x <- c(2.15, 2.95, 1.35, -2.45, 2.45, -0.70, -2.90, -1.35, -2.75)
loadings$label_y <- c(-3.22, 1.75, -0.95, -0.05, 1.05, -3.15, -2.40, 1.75, 2.40)
loadings$label_hjust <- ifelse(loadings$label_x > 0, 0, 1)
loading_scale <- 6.0

p5a <- ggplot(scores, aes(PC1, PC2)) +
  geom_point(aes(fill = entropy, shape = campaign), size = 2.0, colour = "#333333", stroke = 0.25, alpha = 0.88) +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1 * loading_scale, yend = PC2 * loading_scale),
               inherit.aes = FALSE, arrow = arrow(length = unit(1.5, "mm")), colour = "#333333", linewidth = 0.4) +
  geom_text(data = loadings, aes(label_x, label_y, label = label, hjust = label_hjust),
            inherit.aes = FALSE, size = 2.5) +
  scale_fill_gradient(low = "#F7FBFF", high = "#B2182B", limits = c(0, 1), name = "SOM k=3\nentropy") +
  scale_shape_manual(values = c(21, 24, 22, 23), name = "Survey") +
  labs(x = sprintf("PC1 (%.1f%%)", 100 * summary(gradient$full_pca)$importance[2, 1]),
       y = sprintf("PC2 (%.1f%%)", 100 * summary(gradient$full_pca)$importance[2, 2])) + base_theme +
  theme(legend.position = "right")

pca_long <- gradient$pca_stability %>%
  select(procrustes_correlation, pairwise_distance_spearman) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
    procrustes_correlation = "Procrustes correlation",
    pairwise_distance_spearman = "Pairwise-distance correlation"
  ))
p5b <- ggplot(pca_long, aes(metric, value, fill = metric)) +
  geom_violin(alpha = 0.25, colour = NA, width = 0.85) +
  geom_boxplot(width = 0.20, outlier.shape = NA, fill = "white", linewidth = 0.4) +
  geom_point(position = position_jitter(width = 0.09), size = 0.45, alpha = 0.25) +
  scale_fill_manual(values = c("#0072B2", "#009E73")) + coord_cartesian(ylim = c(0.97, 1.002)) +
  labs(x = NULL, y = "Agreement with full-data PCA") + base_theme +
  theme(legend.position = "none", axis.text.x = element_text(angle = 18, hjust = 1))

som_fidelity_plot <- gradient$som_fidelity %>%
  mutate(grid = factor(grid, levels = c("5x5", "7x5", "8x6", "10x7", "12x8")))
p5c <- ggplot(som_fidelity_plot, aes(input_to_map_distance_spearman, component_relation_spearman, colour = grid)) +
  geom_point(size = 1.45, alpha = 0.70) +
  geom_vline(xintercept = median(gradient$som_fidelity$input_to_map_distance_spearman), linetype = "dashed", linewidth = 0.35) +
  geom_hline(yintercept = median(gradient$som_fidelity$component_relation_spearman), linetype = "dashed", linewidth = 0.35) +
  scale_colour_viridis_d(option = "C", end = 0.88) +
  labs(x = "Input-to-map distance correlation", y = "Component-relation correlation", colour = "SOM grid") + base_theme +
  theme(legend.position = "right")

k3_sample <- resampling$sample_metrics %>% filter(k == 3L) %>%
  mutate(method = factor(method, levels = names(method_colours)))
p5d <- ggplot(k3_sample, aes(method, entropy, fill = method)) +
  geom_violin(alpha = 0.23, colour = NA, scale = "width") +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", linewidth = 0.4) +
  geom_point(position = position_jitter(width = 0.08), size = 0.35, alpha = 0.16) +
  scale_fill_manual(values = method_colours) + coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "Normalised assignment entropy (k=3)") + base_theme +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))

fig5 <- (p5a | p5b) / (p5c | p5d) + plot_layout(widths = c(1.15, 0.85)) + plot_annotation(tag_levels = "A")
save_multiformat(fig5, "Fig5_continuous_geometry_and_uncertainty", 180, 158)
write_caption("Fig5_continuous_geometry_and_uncertainty", paste(
  "Persistence of continuous multivariate geometry contrasted with hard-assignment uncertainty.",
  "(A) The full-data PCA plane, with arrows denoting variable loadings and colour denoting SOM k=3 assignment entropy.",
  "(B) Geometry agreement after training-only transformation and PCA in 200 grouped resamples.",
  "(C) Distance and component-relation preservation for 250 SOMs.",
  "(D) Sample assignment entropy at k=3. High PCA agreement does not imply a uniquely recoverable partition."
))

# Fig. 6: campaign transfer and sensitivity ----------------------------------
transfer_plot <- transfer$summary %>% mutate(method = factor(method, levels = names(method_colours)))
p6a <- ggplot(transfer_plot, aes(k, median_ari, colour = held_out_campaign, linetype = held_out_campaign, shape = held_out_campaign, group = held_out_campaign)) +
  geom_hline(yintercept = 0, colour = "#666666", linewidth = 0.3) +
  geom_line(linewidth = 0.65) + geom_point(size = 1.25) +
  facet_wrap(~method, nrow = 2) +
  scale_colour_manual(values = campaign_colours) + scale_linetype_manual(values = campaign_linetypes) +
  scale_shape_manual(values = campaign_shapes) + scale_x_continuous(breaks = K_VALUES) +
  labs(x = "Number of clusters (k)", y = "ARI to full-data consensus", colour = NULL, linetype = NULL, shape = NULL) + base_theme +
  theme(legend.position = "top", legend.text = element_text(size = 6.5), legend.title = element_text(size = 6.5))

p6b <- ggplot(transfer_plot, aes(k, median_distance_ratio, colour = held_out_campaign, linetype = held_out_campaign, shape = held_out_campaign, group = held_out_campaign)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "#555555", linewidth = 0.35) +
  geom_line(linewidth = 0.65) + geom_point(size = 1.25) +
  facet_wrap(~method, nrow = 2) +
  scale_colour_manual(values = campaign_colours) + scale_linetype_manual(values = campaign_linetypes) +
  scale_shape_manual(values = campaign_shapes) + scale_x_continuous(breaks = K_VALUES) +
  labs(x = "Number of clusters (k)", y = "Test / training distance ratio", colour = "Held-out survey", linetype = "Held-out survey", shape = "Held-out survey") + base_theme +
  theme(legend.position = "none")

sens_labels <- c(
  field4 = "Field chemistry", nutrient4 = "Nutrients", nutrient_chla5 = "Chl-a + nutrients",
  raw_core9 = "Untransformed primary", partial5_native_na = "Native missing distance"
)
sens_summary <- sensitivity$module_summary %>%
  left_join(
    readr::read_csv(file.path(TABLE_DIR, "TableS_sensitivity_to_primary_consensus.csv"), show_col_types = FALSE) %>%
      mutate(analysis = recode(analysis, partial5 = "partial5_native_na")),
    by = c("analysis", "k")
  ) %>%
  mutate(analysis_label = factor(sens_labels[analysis], levels = unname(sens_labels))) %>%
  select(analysis_label, k, median_ari_to_consensus, ari_to_primary_som_consensus) %>%
  pivot_longer(c(median_ari_to_consensus, ari_to_primary_som_consensus), names_to = "metric", values_to = "ari") %>%
  mutate(metric = recode(metric,
    median_ari_to_consensus = "Within-variant stability",
    ari_to_primary_som_consensus = "Agreement with primary"
  ))
p6c <- ggplot(sens_summary, aes(k, ari, colour = analysis_label, linetype = metric, shape = metric, group = interaction(analysis_label, metric))) +
  geom_line(linewidth = 0.65) + geom_point(size = 1.25) +
  scale_x_continuous(breaks = K_VALUES) + scale_colour_viridis_d(option = "D", end = 0.90) +
  scale_linetype_manual(values = c(`Within-variant stability` = "solid", `Agreement with primary` = "dashed")) +
  scale_shape_manual(values = c(`Within-variant stability` = 16, `Agreement with primary` = 1)) +
  labs(x = "Number of clusters (k)", y = "Adjusted Rand index", colour = "Sensitivity analysis", linetype = NULL, shape = NULL) + base_theme +
  theme(legend.position = "right", legend.text = element_text(size = 6.5))

exclusion_plot <- sensitivity$exclusion %>%
  mutate(method = factor(method, levels = rev(names(method_colours))))
p6d <- ggplot(exclusion_plot, aes(k, method, fill = ari_to_primary_consensus)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f", ari_to_primary_consensus)), size = 2.25) +
  scale_fill_gradient2(low = "#B2182B", mid = "#F7F7F7", high = "#2166AC", midpoint = 0.5, limits = c(0, 1),
                       name = "ARI") +
  scale_x_continuous(breaks = K_VALUES) +
  labs(x = "Number of clusters (k)", y = NULL, title = "Exclusion of the 2024-2025 survey") + base_theme +
  theme(panel.grid = element_blank(), legend.position = "right")

fig6 <- (p6a | p6b) / (p6c | p6d) + plot_annotation(tag_levels = "A")
save_multiformat(fig6, "Fig6_transfer_and_sensitivity", 180, 168)
write_caption("Fig6_transfer_and_sensitivity", paste(
  "Transfer and sensitivity analyses.",
  "(A) Agreement between leave-one-survey-out predictions and the corresponding full-data consensus labels; these labels are a concordance reference, not ground truth.",
  "(B) Median held-out-to-training centroid-distance ratio as a domain-shift diagnostic.",
  "(C) Within-variant SOM ensemble stability and agreement with the primary consensus for variable modules, untransformed inputs and native missing-distance analysis.",
  "(D) Agreement after excluding the 2024–2025 survey."
))

# Supplementary figures -------------------------------------------------------
consensus_frames <- lapply(names(method_colours), function(method) {
  cons <- resampling$consensus[[method]]$k3
  order_index <- order(cons$labels)
  mat <- cons$coassociation[order_index, order_index]
  data.frame(
    method = method,
    x = rep(seq_len(nrow(mat)), each = nrow(mat)),
    y = rep(seq_len(nrow(mat)), times = nrow(mat)),
    coassociation = as.vector(mat)
  )
}) %>% bind_rows()
s1 <- ggplot(consensus_frames, aes(x, y, fill = coassociation)) +
  geom_raster() + facet_wrap(~method, nrow = 1) +
  scale_fill_gradient(low = "white", high = "#2166AC", limits = c(0, 1)) +
  coord_equal() + labs(x = "Samples ordered by consensus cluster", y = NULL, fill = "Coassociation") + base_theme +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
save_multiformat(s1, "FigS1_k3_consensus_matrices", 180, 62, FIG_SUPP_DIR)
write_caption("FigS1_k3_consensus_matrices", "Coassociation matrices for the k=3 resampling consensus of each model family. Samples are ordered within each panel by its own consensus labels.", FIG_SUPP_DIR)

gmm_count <- resampling$gmm_models %>% dplyr::count(k, model_name)
s2 <- ggplot(gmm_count, aes(k, model_name, fill = n)) +
  geom_tile(colour = "white") + geom_text(aes(label = n), size = 2.4) +
  scale_fill_gradient(low = "#F7FBFF", high = "#2166AC") + scale_x_continuous(breaks = K_VALUES) +
  labs(x = "Number of clusters (k)", y = "GMM covariance model", fill = "Resamples") + base_theme
save_multiformat(s2, "FigS2_gmm_covariance_selection", 120, 75, FIG_SUPP_DIR)
write_caption("FigS2_gmm_covariance_selection", "GMM covariance structures selected by Bayesian information criterion within each of 200 training resamples.", FIG_SUPP_DIR)

som_agreement_with_grid <- ensemble$model_agreement %>%
  left_join(ensemble$metrics %>% select(model, grid), by = "model")
s3 <- ggplot(som_agreement_with_grid, aes(factor(k), ari_to_consensus, fill = grid)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.35) +
  scale_fill_viridis_d(option = "C", end = 0.88) +
  labs(x = "Number of clusters (k)", y = "ARI to 250-model SOM consensus", fill = "SOM grid") + base_theme +
  theme(legend.position = "bottom")
save_multiformat(s3, "FigS3_som_grid_seed_partition_agreement", 150, 92, FIG_SUPP_DIR)
write_caption("FigS3_som_grid_seed_partition_agreement", "Agreement of every SOM grid-seed partition with the corresponding 250-model consensus.", FIG_SUPP_DIR)

validity_long <- resampling$internal_validity %>%
  pivot_longer(c(silhouette, calinski_harabasz, davies_bouldin), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, silhouette = "Silhouette", calinski_harabasz = "Calinski-Harabasz", davies_bouldin = "Davies-Bouldin"))
s4 <- ggplot(validity_long, aes(k, value, colour = method, linetype = method, shape = method)) +
  geom_line(linewidth = 0.65) + geom_point(size = 1.25) + facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = method_colours) +
  scale_linetype_manual(values = c(SOM = "solid", `K-means` = "dashed", `Ward.D2` = "dotdash", GMM = "dotted")) +
  scale_shape_manual(values = c(SOM = 16, `K-means` = 17, `Ward.D2` = 15, GMM = 18)) +
  scale_x_continuous(breaks = K_VALUES) +
  labs(x = "Number of clusters (k)", y = NULL, colour = NULL, linetype = NULL, shape = NULL) + base_theme + theme(legend.position = "bottom")
save_multiformat(s4, "FigS4_internal_validity_indices", 180, 72, FIG_SUPP_DIR)
write_caption("FigS4_internal_validity_indices", "Internal geometric indices for each full resampling-consensus partition. These indices are auxiliary and are not treated as ecological truth or stability thresholds.", FIG_SUPP_DIR)

loading_summary <- gradient$loading_stability %>%
  group_by(variable) %>% summarise(median = median(loading_vector_cosine), lower = quantile(loading_vector_cosine, 0.025), upper = quantile(loading_vector_cosine, 0.975), .groups = "drop")
s5 <- ggplot(loading_summary, aes(reorder(variable, median), median)) +
  geom_hline(yintercept = 0, colour = "#777777", linewidth = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.18, linewidth = 0.45) +
  geom_point(size = 2.0, colour = "#0072B2") + coord_flip() +
  labs(x = NULL, y = "PCA loading-vector cosine after alignment") + base_theme
save_multiformat(s5, "FigS5_pca_loading_stability", 120, 92, FIG_SUPP_DIR)
write_caption("FigS5_pca_loading_stability", "Stability of two-dimensional PCA loading vectors across 200 grouped resamples after orthogonal Procrustes alignment.", FIG_SUPP_DIR)


message("Generated Figs. 3-6 and supplementary figures. Figs. 1-2 are outside this public package.")
