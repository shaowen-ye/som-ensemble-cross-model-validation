# Recreate manuscript Figs. 3-6 and the analytical supplementary figures S1-S7.
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
# First-round revision outputs (scripts 10-14).
topology <- readRDS(file.path(MODEL_DIR, "r1_som_topology.rds"))
calibration <- readRDS(file.path(MODEL_DIR, "r1_transfer_calibration.rds"))
like_for_like <- readRDS(file.path(MODEL_DIR, "r1_sensitivity_like_for_like.rds"))

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
  geom_point(position = position_jitter(width = 0.10, height = 0, seed = 3), size = 0.45, alpha = 0.35) +
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
  guides(colour = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
  labs(x = "Number of clusters (k)", y = "Normalised assignment entropy", colour = NULL, fill = NULL) + base_theme +
  theme(legend.position = "inside", legend.position.inside = c(0.5, 0.995), legend.justification = c(0.5, 1),
        legend.background = element_rect(fill = "white", colour = NA), legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6.5))

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
  geom_point(position = position_jitter(width = 0.09, seed = 5), size = 0.45, alpha = 0.25) +
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
# The median printed above each box is the value quoted in the text; the
# diamond marks the mean, which the right-skewed distributions pull upwards.
k3_centre <- k3_sample %>% group_by(method) %>%
  summarise(median = median(entropy), mean = mean(entropy), .groups = "drop")
p5d <- ggplot(k3_sample, aes(method, entropy, fill = method)) +
  geom_violin(alpha = 0.23, colour = NA, scale = "width") +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", linewidth = 0.4) +
  geom_point(position = position_jitter(width = 0.08, seed = 7), size = 0.35, alpha = 0.16) +
  geom_point(data = k3_centre, aes(method, mean), inherit.aes = FALSE,
             shape = 23, size = 1.7, fill = "white", colour = "#222222", stroke = 0.45) +
  geom_text(data = k3_centre, aes(method, 1.075, label = sprintf("median\n%.3f", median)),
            inherit.aes = FALSE, size = 2.15, lineheight = 0.85) +
  scale_fill_manual(values = method_colours) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1.13)) +
  labs(x = NULL, y = "Normalised assignment entropy (k=3)") + base_theme +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))

fig5 <- (p5a | p5b) / (p5c | p5d) + plot_layout(widths = c(1.15, 0.85)) + plot_annotation(tag_levels = "A")
save_multiformat(fig5, "Fig5_continuous_geometry_and_uncertainty", 180, 158)
write_caption("Fig5_continuous_geometry_and_uncertainty", paste(
  "Persistence of continuous multivariate geometry contrasted with hard-assignment uncertainty.",
  "(A) The full-data PCA plane, with arrows denoting variable loadings and colour denoting SOM k=3 assignment entropy.",
  "(B) Geometry agreement after training-only transformation and PCA in 200 grouped resamples.",
  "(C) Distance and component-relation preservation for 250 SOMs.",
  "(D) Sample assignment entropy at k=3; printed values are medians and diamonds are means. High PCA agreement does not imply a uniquely recoverable partition."
))

# Fig. 6: survey transfer and sensitivity, judged against matched references ---
# Every observed value is shown beside the range produced by removing the same
# number of randomly chosen events (scripts 13 and 14), so that sample-size
# effects are not read as survey-domain effects.
method_dodge <- position_dodge(width = 0.74)
transfer_k3 <- calibration$comparison %>%
  filter(k == 3L, held_out_campaign != "2024-2025") %>%
  mutate(method = factor(method, levels = names(method_colours)),
         survey = factor(sprintf("%s\n(n = %d)", held_out_campaign, n_test)))
matched_layers <- function(lower, upper, middle, observed) {
  list(
    geom_linerange(aes(ymin = {{ lower }}, ymax = {{ upper }}, colour = method), position = method_dodge,
                   linewidth = 2.3, alpha = 0.28),
    geom_errorbar(aes(ymin = {{ middle }}, ymax = {{ middle }}, colour = method), position = method_dodge,
                  width = 0.17, linewidth = 0.45),
    geom_point(aes(y = {{ observed }}, fill = method), position = method_dodge, shape = 21, size = 1.9,
               colour = "#222222", stroke = 0.3)
  )
}
p6a <- ggplot(transfer_k3, aes(survey, group = method)) +
  geom_hline(yintercept = 0, colour = "#777777", linewidth = 0.3) +
  matched_layers(ari_lower, ari_upper, ari_median, loso_ari) +
  scale_colour_manual(values = method_colours) + scale_fill_manual(values = method_colours) +
  coord_cartesian(ylim = c(-0.05, 1.02)) +
  labs(x = "Held-out survey", y = "ARI with full-data consensus (k = 3)", colour = NULL, fill = NULL) + base_theme +
  theme(legend.position = "top", legend.text = element_text(size = 6.5))

p6b <- ggplot(transfer_k3, aes(survey, group = method)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "#555555", linewidth = 0.35) +
  matched_layers(ratio_lower, ratio_upper, ratio_median, loso_ratio) +
  scale_colour_manual(values = method_colours) + scale_fill_manual(values = method_colours) +
  labs(x = "Held-out survey", y = "Test / training distance ratio (k = 3)", colour = NULL, fill = NULL) + base_theme +
  theme(legend.position = "none")

module_levels <- c("Field chemistry", "Nutrients", "Chl-a + nutrients", "Untransformed nine variables", "Native missing distance")
module_k3 <- like_for_like$modules %>%
  filter(k == 3L, analysis %in% module_levels) %>%
  mutate(events = ifelse(cohort == "Primary 390 events", "Same 390 events", "Larger module cohort"),
         analysis = factor(recode(analysis, `Untransformed nine variables` = "Untransformed",
                                  `Native missing distance` = "Missing distance"),
                           levels = rev(c("Field chemistry", "Nutrients", "Chl-a + nutrients", "Untransformed", "Missing distance")))) %>%
  select(analysis, events, within_median, agreement_primary_grid_seed) %>%
  pivot_longer(c(within_median, agreement_primary_grid_seed), names_to = "metric", values_to = "ari") %>%
  mutate(metric = factor(recode(metric, within_median = "Within-variant stability",
                                agreement_primary_grid_seed = "Agreement with primary"),
                         levels = c("Within-variant stability", "Agreement with primary")))
primary_reference <- data.frame(
  metric = factor("Within-variant stability", levels = levels(module_k3$metric)),
  value = like_for_like$baseline$within_median[like_for_like$baseline$k == 3L]
)
p6c <- ggplot(module_k3, aes(ari, analysis, shape = events, fill = events)) +
  geom_vline(data = primary_reference, aes(xintercept = value), linetype = "dashed", colour = "#555555", linewidth = 0.35) +
  geom_point(size = 2.0, colour = "#333333", stroke = 0.35) +
  facet_wrap(~metric, nrow = 1) +
  scale_shape_manual(values = c(`Same 390 events` = 21, `Larger module cohort` = 24)) +
  scale_fill_manual(values = c(`Same 390 events` = "#0072B2", `Larger module cohort` = "white")) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = "Adjusted Rand index (k = 3)", y = NULL, shape = NULL, fill = NULL) + base_theme +
  theme(legend.position = "top", legend.text = element_text(size = 6.5))

exclusion_plot <- like_for_like$exclusion %>% mutate(method = factor(method, levels = names(method_colours)))
p6d <- ggplot(exclusion_plot, aes(k)) +
  geom_ribbon(aes(ymin = random_omission_lower, ymax = random_omission_upper, fill = method), alpha = 0.18) +
  geom_line(aes(y = random_omission_median, colour = method), linetype = "dashed", linewidth = 0.35) +
  geom_linerange(data = exclusion_plot %>% filter(method == "SOM"),
                 aes(x = k + 0.22, ymin = single_map_lower, ymax = single_map_upper), colour = "#8C8C8C", linewidth = 0.55) +
  geom_line(aes(y = ensemble_consensus_ari, colour = method), linewidth = 0.7) +
  geom_point(aes(y = ensemble_consensus_ari, colour = method), size = 1.2) +
  facet_wrap(~method, nrow = 2) +
  scale_colour_manual(values = method_colours) + scale_fill_manual(values = method_colours) +
  scale_x_continuous(breaks = K_VALUES) + coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Number of clusters (k)", y = "ARI with primary consensus,\n2024-2025 survey excluded") + base_theme +
  theme(legend.position = "none")

fig6 <- (p6a | p6b) / (p6c | p6d) + plot_annotation(tag_levels = "A")
save_multiformat(fig6, "Fig6_transfer_and_sensitivity", 180, 168)
write_caption("Fig6_transfer_and_sensitivity", paste(
  "Survey transfer and sensitivity analyses judged against matched references.",
  "(A) Agreement at k = 3 between leave-one-survey-out predictions and the full-data consensus (circles), beside the 95% range (bars) and median (ticks) of 50 random omissions of the same number of events. The 2024–2025 survey (16 events) is not shown: at k = 3 its SOM, K-means and Ward.D2 references are single clusters, and its GMM reference splits the events 12/4 (Table S20). Consensus labels are a concordance reference, not ground truth.",
  "(B) Held-out-to-training centroid-distance ratio, shown in the same way; the dashed line marks a ratio of one.",
  "(C) Variable modules and preprocessing at k = 3: within-variant stability of each 50-map SOM ensemble and its agreement with the primary grid–seed consensus (the 50 primary maps with the same grids and seeds), fitted to the same 390 events or to the larger module cohort; the dashed line marks the stability of the primary nine-variable ensemble.",
  "(D) Agreement with the primary consensus after excluding the 2024–2025 survey (solid lines; SOM: 50-map consensus), with the median (dashed lines) and 95% range (shading) of 50 random omissions of 16 events; grey bars show the 95% range of single SOM maps.",
  "All-k results are given in Tables S11, S20 and S21 and Fig. S7."
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

# Fig. S6: reproducibility of the SOM map layout (script 10) -------------------
grid_levels <- c("5x5", "7x5", "8x6", "10x7", "12x8")
concordance_plot <- bind_rows(
  topology$full_pairs %>% filter(grid_1 == grid_2) %>%
    transmute(group = grid_1, value = map_distance_spearman, fits = "Full-data SOMs"),
  topology$resample_pairs %>% filter(grid_1 == grid_2) %>%
    transmute(group = "Resample\nSOMs", value = map_distance_spearman, fits = "Resample SOMs"),
  topology$quantised_pairs %>% filter(grid_1 == grid_2) %>%
    transmute(group = "PCA on\nSOM grid", value = map_distance_spearman, fits = "Resample PCA on SOM grid"),
  topology$pca_pairs %>%
    transmute(group = "Resample\nPCA", value = map_distance_spearman, fits = "Resample PCA")
) %>% mutate(group = factor(group, levels = c(grid_levels, "Resample\nSOMs", "PCA on\nSOM grid", "Resample\nPCA")),
             fits = factor(fits, levels = c("Full-data SOMs", "Resample SOMs", "Resample PCA on SOM grid", "Resample PCA")))
fit_colours <- c(`Full-data SOMs` = "#0072B2", `Resample SOMs` = "#56B4E9",
                 `Resample PCA on SOM grid` = "#E69F00", `Resample PCA` = "#009E73")
s6a <- ggplot(concordance_plot, aes(group, value, fill = fits)) +
  geom_violin(scale = "width", alpha = 0.35, colour = NA) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", linewidth = 0.35) +
  scale_fill_manual(values = fit_colours) + coord_cartesian(ylim = c(0, 1)) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = "Same-grid pairs of fits", y = "Map-layout concordance (Spearman)", fill = NULL) + base_theme +
  theme(legend.position = "top", legend.text = element_text(size = 6.5))
overlap_plot <- bind_rows(
  topology$full_overlap %>% transmute(group = grid, observed = neighbour_overlap, chance = chance_overlap),
  topology$resample_overlap %>% transmute(group = "Resample\nSOMs", observed = neighbour_overlap, chance = chance_overlap)
) %>%
  pivot_longer(c(observed, chance), names_to = "type", values_to = "value") %>%
  mutate(group = factor(group, levels = c(grid_levels, "Resample\nSOMs")),
         type = factor(recode(type, observed = "Observed", chance = "Events relabelled at random"), levels = c("Observed", "Events relabelled at random")))
s6b <- ggplot(overlap_plot, aes(group, value, fill = type)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.35, width = 0.6, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c(Observed = "#0072B2", `Events relabelled at random` = "#D9D9D9")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Same-grid pairs of fits", y = "Map-neighbour overlap (Jaccard)", fill = NULL) + base_theme +
  theme(legend.position = "top", legend.text = element_text(size = 6.5))
s6 <- (s6a | s6b) + plot_layout(widths = c(1.15, 1)) + plot_annotation(tag_levels = "A")
save_multiformat(s6, "FigS6_som_map_layout_reproducibility", 180, 82, FIG_SUPP_DIR)
write_caption("FigS6_som_map_layout_reproducibility", "Reproducibility of the SOM map layout. (A) Map-layout concordance: Spearman correlation between the grid distances separating the best matching units of all event pairs in two SOM fits, for same-grid pairs of the 250 full-data maps and of the 200 grouped-resample maps. For reference, the aligned two-dimensional PCA scores of the same resamples were quantised onto the matching SOM grid and compared in the same way, and event distances in the continuous PCA plane are shown last. (B) Mean event-level Jaccard overlap of map neighbours (events mapped to the same or an adjacent unit) between two fits, beside the overlap expected after relabelling events at random.", FIG_SUPP_DIR)

# Fig. S7: survey transfer against matched-size random omission, all k ---------
transfer_all <- calibration$comparison %>%
  mutate(method = factor(method, levels = names(method_colours)),
         held_out = sprintf("%s (n = %d)", held_out_campaign, n_test),
         reference = factor(ifelse(reference_clusters > 1L, "multi", "single"), levels = c("multi", "single")))
s7 <- ggplot(transfer_all, aes(k)) +
  geom_ribbon(aes(ymin = ari_lower, ymax = ari_upper), fill = "#BDBDBD", alpha = 0.55) +
  geom_line(aes(y = ari_median), linetype = "dashed", colour = "#555555", linewidth = 0.35) +
  geom_line(aes(y = loso_ari, colour = method), linewidth = 0.65) +
  geom_point(aes(y = loso_ari, colour = method, shape = reference), size = 1.1, fill = "white", stroke = 0.45) +
  facet_grid(method ~ held_out) +
  scale_colour_manual(values = method_colours) + scale_shape_manual(values = c(multi = 16, single = 21)) +
  scale_x_continuous(breaks = K_VALUES) +
  coord_cartesian(ylim = c(-0.1, 1)) +
  labs(x = "Number of clusters (k)", y = "ARI with full-data consensus") + base_theme +
  theme(legend.position = "none", strip.text = element_text(size = 6.5))
save_multiformat(s7, "FigS7_transfer_matched_omission_all_k", 180, 130, FIG_SUPP_DIR)
write_caption("FigS7_transfer_matched_omission_all_k", "Leave-one-survey-out agreement (coloured lines) compared with 50 random omissions of the same number of events (grey: 95% range; dashed: median) for every method and k. Open symbols mark combinations in which the method's own reference places all held-out events in one cluster, so that ARI can only be 0 or 1. The 2024–2025 survey contributes 16 events and is not interpreted.", FIG_SUPP_DIR)

message("Generated Figs. 3-6 and supplementary figures. Figs. 1-2 are outside this public package.")
