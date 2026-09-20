# ==============================================================================
# 09_figures.R
# Purpose: produce publication figures for the manuscript.
# All figures are built from intermediate objects produced by Stages 1-4.
# Figures saved to manuscript/figures/ as PDF (vector) + PNG (raster preview).
# ==============================================================================

source(here::here("R", "00_setup.R"))
STAGE <- "figures"

ad_bulk <- load_intermediate("ad_bulk")
gbm_bulk <- load_intermediate("gbm_bulk")
stage1 <- load_intermediate("stage1_ad_vulnerable")
stage2 <- load_intermediate("stage2_gbm_neural_mimicry")
stage3 <- load_intermediate("stage3_convergence")
stage4 <- load_intermediate("stage4_bulk_reanalysis")

fig_dir <- here::here(PARAMS$paths$figures)

save_figure <- function(plot, name, w = 6, h = 4) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), plot, width = w, height = h)
  ggsave(file.path(fig_dir, paste0(name, ".png")), plot, width = w, height = h, dpi = 300)
  log_msg(STAGE, sprintf("Saved: %s", name))
}

# Theme
theme_pub <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", color = NA),
        legend.position = "right")

# ------------------------------------------------------------------------------
# Figure 1. Cell composition shifts (BRETIGEA scores)
# ------------------------------------------------------------------------------

ad_scores <- stage4$ad_bretigea %>%
  tibble::rownames_to_column("sample_id") %>%
  left_join(ad_bulk$meta %>% select(sample_id, diagnosis), by = "sample_id") %>%
  pivot_longer(cols = -c(sample_id, diagnosis), names_to = "cell_type", values_to = "score") %>%
  mutate(cohort = "AD")

gbm_scores <- stage4$gbm_bretigea %>%
  tibble::rownames_to_column("sample_id") %>%
  left_join(gbm_bulk$meta %>% select(sample_id, diagnosis), by = "sample_id") %>%
  pivot_longer(cols = -c(sample_id, diagnosis), names_to = "cell_type", values_to = "score") %>%
  mutate(cohort = "GBM")

fig1_data <- bind_rows(ad_scores, gbm_scores) %>%
  mutate(cell_type = factor(cell_type,
                             levels = c("neu", "ast", "mic", "oli", "opc", "end"),
                             labels = c("Neurons", "Astrocytes", "Microglia",
                                        "Oligodendrocytes", "OPCs", "Endothelial")),
         group = ifelse(diagnosis %in% c("AD", "GBM"), "Disease", "Control"))

p_fig1 <- ggplot(fig1_data, aes(x = group, y = score, fill = group)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  facet_grid(cohort ~ cell_type, scales = "free_y") +
  scale_fill_manual(values = c("Control" = "#4A90D9", "Disease" = "#D94A4A")) +
  labs(x = NULL, y = "BRETIGEA cell-type score") +
  theme_pub + theme(legend.position = "none")

save_figure(p_fig1, "fig1_cell_composition", w = 8, h = 5)

# ------------------------------------------------------------------------------
# Figure 2. DEG collapse under cell composition control (H1)
# ------------------------------------------------------------------------------

deg_data <- stage4$shared_counts %>%
  mutate(model_label = c("A" = "Uncorrected",
                         "B" = "BRETIGEA-controlled",
                         "C" = "bMIND-controlled")[model]) %>%
  filter(!is.na(n_shared))

p_fig2 <- ggplot(deg_data, aes(x = model_label, y = n_shared)) +
  geom_col(fill = "#4A90D9", width = 0.6) +
  geom_text(aes(label = n_shared), vjust = -0.4, size = 3.5) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0, 10, 100, 500, 1000, 5000),
                     expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Shared AD-GBM DEGs (log scale)",
       title = sprintf("H1 decision: %s", stage4$h1_decision)) +
  theme_pub

save_figure(p_fig2, "fig2_h1_deg_collapse", w = 5, h = 4)

# ------------------------------------------------------------------------------
# Figure 3. RRHO2 heatmap (Test 1 of H2)
# ------------------------------------------------------------------------------

# The RRHO2 object contains the hypermat; visualise it as a heatmap
if (!is.null(stage3$rrho_object$hypermat)) {
  hm <- stage3$rrho_object$hypermat
  hm_df <- as.data.frame(as.table(hm)) %>%
    rename(x = Var1, y = Var2, value = Freq) %>%
    mutate(x = as.integer(x), y = as.integer(y))

  p_fig3 <- ggplot(hm_df, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    scale_fill_gradient2(low = "white", mid = "#FFDDBB", high = "#B30000",
                         midpoint = median(hm_df$value, na.rm = TRUE),
                         name = "-log10(P)") +
    labs(x = "AD signed -log10(P) rank",
         y = "GBM signed -log10(P) rank",
         title = sprintf("RRHO2: max -log10(P) = %.2f (BH P = %.3g)",
                         stage3$rrho_test$max_neg_log10_p,
                         stage3$rrho_test$max_p_bh)) +
    theme_pub

  save_figure(p_fig3, "fig3_rrho2", w = 6, h = 5)
}

# ------------------------------------------------------------------------------
# Figure 4. Permutation null distribution (Test 3 of H2)
# ------------------------------------------------------------------------------

perm_df <- data.frame(odds_ratio = stage3$perm_or_distribution)
p_fig4 <- ggplot(perm_df, aes(x = odds_ratio)) +
  geom_histogram(bins = 50, fill = "#CCCCCC", color = "black", size = 0.2) +
  geom_vline(xintercept = stage3$permutation_test$observed_or,
             color = "#D94A4A", linetype = "dashed", size = 1) +
  annotate("text", x = stage3$permutation_test$observed_or,
           y = Inf, vjust = 1.5, hjust = -0.1,
           label = sprintf("Observed OR = %.2f\nP = %.4f",
                           stage3$permutation_test$observed_or,
                           stage3$permutation_test$p_empirical),
           color = "#D94A4A", size = 3) +
  labs(x = "Odds ratio (random gene sets)",
       y = "Frequency",
       title = "Permutation test: null distribution vs observed") +
  theme_pub

save_figure(p_fig4, "fig4_permutation_null", w = 5, h = 4)

# ------------------------------------------------------------------------------
# Figure 5. Signature overlap Venn / summary
# ------------------------------------------------------------------------------

# Simple bar summary of overlap
overlap_summary <- data.frame(
  category = c("AD-vulnerable\nonly",
               sprintf("Overlap\n(n=%d)", stage3$hypergeometric_test$overlap_n),
               "GBM-mimicry\nonly"),
  count = c(PARAMS$stage3$top_n_hypergeometric - stage3$hypergeometric_test$overlap_n,
            stage3$hypergeometric_test$overlap_n,
            PARAMS$stage3$top_n_hypergeometric - stage3$hypergeometric_test$overlap_n)
)

p_fig5 <- ggplot(overlap_summary, aes(x = category, y = count, fill = category)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = count), vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = c("#4A90D9", "#D94A4A", "#7CB342")) +
  labs(x = NULL, y = sprintf("Top %d genes", PARAMS$stage3$top_n_hypergeometric)) +
  theme_pub + theme(legend.position = "none")

save_figure(p_fig5, "fig5_signature_overlap", w = 5, h = 4)

snapshot_session(STAGE)
log_msg(STAGE, "Figures complete.")
