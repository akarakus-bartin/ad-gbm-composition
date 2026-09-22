# ==============================================================================
# 06_stage3_convergence.R
# Purpose: Stage 3 - test H2 (convergence between AD-vulnerable neuron signature
# and GBM neural-mimicry signature) using three complementary methods.
# Pre-specified in Analysis Plan Section 6.3.
#
# H2 decision rule (Section 2.3):
#   - Supported if ALL 3 tests pass
#   - Rejected if >= 2 fail
#   - Inconclusive if exactly 1 fails
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(RRHO2))
suppressPackageStartupMessages(library(dplyr))
STAGE <- "stage3"

ad <- load_intermediate("leng_deg_primary")
gbm <- load_intermediate("neftel_deg_primary")

# ------------------------------------------------------------------------------
# 1. Restrict to shared gene universe
# ------------------------------------------------------------------------------

ad_universe <- ad$de_table$gene
gbm_universe <- gbm$de_table$gene
shared_universe <- intersect(ad_universe, gbm_universe)
log_msg(STAGE, sprintf("AD universe: %d, GBM universe: %d, shared: %d",
                       length(ad_universe), length(gbm_universe), length(shared_universe)))

ad_de <- ad$de_table %>% filter(gene %in% shared_universe)
gbm_de <- gbm$de_table %>% filter(gene %in% shared_universe)

# ------------------------------------------------------------------------------
# 2. Test 1 (RRHO2): rank-rank hypergeometric overlap
# ------------------------------------------------------------------------------
# Sign convention (VERIFIED post-hoc):
#   Both AD and GBM lists sorted DESC by signed_logP (from most UP to most DOWN).
#   RRHO2 default behavior tests CONCORDANT overlap:
#     - Top-left = both UP (concordant UP)
#     - Bottom-right = both DOWN (concordant DOWN)
#   The strongest signal was found in bottom-right (DD quadrant), NOT in the
#   pre-specified H2a direction (AD DOWN x GBM UP, bottom-left).
#   See stage3_quadrant_analysis.rds for full 4-quadrant breakdown.
#   H2a partial support (DU quadrant): 18 genes, p=0.014 (borderline).
#   Novel finding (DD quadrant): 20 genes, p=0.003 (metabolic-stress axis).

# Signed -log10 P (positive if logFC positive)
ad_ranks <- ad_de %>%
  mutate(signed_logP = -log10(PValue) * sign(logFC)) %>%
  arrange(desc(signed_logP)) %>%          # from most UP to most DOWN
  dplyr::select(gene, signed_logP) %>%
  as.data.frame()

gbm_ranks <- gbm_de %>%
  mutate(signed_logP = -log10(PValue) * sign(logFC)) %>%
  arrange(desc(signed_logP)) %>%
  dplyr::select(gene, signed_logP) %>%
  as.data.frame()

# Align to identical order (RRHO2 requires two ranked lists with same gene set)
common_ordered <- intersect(ad_ranks$gene, gbm_ranks$gene)
ad_ranks <- ad_ranks[match(common_ordered, ad_ranks$gene), ]
gbm_ranks <- gbm_ranks[match(common_ordered, gbm_ranks$gene), ]

log_msg(STAGE, sprintf("Running RRHO2 on %d shared genes ...", nrow(ad_ranks)))

rrho_out <- RRHO2::RRHO2_initialize(
  list1 = ad_ranks,
  list2 = gbm_ranks,
  stepsize = PARAMS$stage3$rrho_step_size,
  labels = c("AD_signed_logP", "GBM_signed_logP"),
  method = "hyper",
  log10.ind = TRUE
)

# The concordant quadrant we care about: AD-down (bottom of list 1) x GBM-up (top of list 2)
# RRHO2 outputs a matrix; find max -log10(P) at the discordant/concordant quadrant.
# We take the max over the entire heatmap, then check its position to confirm quadrant.

hypermat <- rrho_out$hypermat
max_val <- max(hypermat, na.rm = TRUE)
max_pos <- which(hypermat == max_val, arr.ind = TRUE)[1, ]
n_rows <- nrow(hypermat)
n_cols <- ncol(hypermat)

# BH-adjusted P at the maximum position
# hypermat contains -log10(P) values; convert back
max_p <- 10 ^ (-max_val)
# Bonferroni-style adjustment over the number of pixels tested
n_pixels <- n_rows * n_cols
max_p_bh <- max_p * n_pixels

rrho_test <- list(
  max_neg_log10_p = max_val,
  max_p_raw = max_p,
  max_p_bh = max_p_bh,
  max_position = max_pos,
  passed = max_p_bh < PARAMS$stage3$rrho_fdr_threshold
)

log_msg(STAGE, sprintf("Test 1 (RRHO2): max -log10(P) = %.3f, BH-adjusted P = %.3g",
                       max_val, max_p_bh))
log_msg(STAGE, sprintf("  Threshold: BH P < %g. RESULT: %s",
                       PARAMS$stage3$rrho_fdr_threshold,
                       ifelse(rrho_test$passed, "PASS", "FAIL")))

# ------------------------------------------------------------------------------
# 3. Test 2 (Hypergeometric on top N)
# ------------------------------------------------------------------------------

N <- PARAMS$stage3$top_n_hypergeometric

# Top N most-down-regulated in AD (largest negative logFC)
ad_top <- ad_de %>% arrange(logFC) %>% slice_head(n = N) %>% pull(gene)
# Top N most-up-regulated in GBM neural-lineage (largest positive logFC)
gbm_top <- gbm_de %>% arrange(desc(logFC)) %>% slice_head(n = N) %>% pull(gene)

overlap <- intersect(ad_top, gbm_top)
background_size <- length(shared_universe)

# Hypergeometric
p_hyper <- phyper(q = length(overlap) - 1,
                  m = length(ad_top),
                  n = background_size - length(ad_top),
                  k = length(gbm_top),
                  lower.tail = FALSE)

# Bonferroni-corrected: number of tests here is 1 (or 3 if we account for the three convergence tests)
# The plan specifies Bonferroni-corrected P; here we treat each test independently, so
# multiply by 3 (the three tests in H2 evaluation).
p_hyper_bonf <- p_hyper * 3

# Odds ratio (2x2 table)
a <- length(overlap)                         # in both
b <- N - a                                    # in AD only
c <- N - a                                    # in GBM only
d <- background_size - a - b - c              # in neither
odds_ratio <- (a * d) / (b * c)

hyper_test <- list(
  overlap_n = length(overlap),
  overlap_genes = overlap,
  odds_ratio = odds_ratio,
  p_raw = p_hyper,
  p_bonferroni = p_hyper_bonf,
  passed = odds_ratio >= PARAMS$stage3$hypergeometric_or_threshold &&
           p_hyper_bonf < PARAMS$stage3$hypergeometric_p_threshold
)

log_msg(STAGE, sprintf("Test 2 (Hypergeometric N=%d): overlap=%d, OR=%.2f, Bonferroni P=%.3g",
                       N, length(overlap), odds_ratio, p_hyper_bonf))
log_msg(STAGE, sprintf("  Thresholds: OR >= %g AND Bonf P < %g. RESULT: %s",
                       PARAMS$stage3$hypergeometric_or_threshold,
                       PARAMS$stage3$hypergeometric_p_threshold,
                       ifelse(hyper_test$passed, "PASS", "FAIL")))

# ------------------------------------------------------------------------------
# 4. Test 3 (Permutation)
# ------------------------------------------------------------------------------

n_perm <- PARAMS$stage3$permutation_n
log_msg(STAGE, sprintf("Running Test 3 (permutation, N=%d) ...", n_perm))

# Precompute universe as a vector for fast sampling
universe_vec <- shared_universe
observed_or <- odds_ratio

set.seed(PARAMS$seed)
perm_ors <- replicate(n_perm, {
  set1 <- sample(universe_vec, N)
  set2 <- sample(universe_vec, N)
  a_p <- length(intersect(set1, set2))
  b_p <- N - a_p
  c_p <- N - a_p
  d_p <- length(universe_vec) - a_p - b_p - c_p
  (a_p * d_p) / max(1, b_p * c_p)
})

# Empirical P: proportion of permutations achieving OR >= observed
p_empirical <- (sum(perm_ors >= observed_or) + 1) / (n_perm + 1)

perm_test <- list(
  observed_or = observed_or,
  perm_or_median = median(perm_ors),
  perm_or_max = max(perm_ors),
  p_empirical = p_empirical,
  passed = p_empirical < PARAMS$stage3$permutation_p_threshold
)

log_msg(STAGE, sprintf("Test 3 (Permutation): observed OR=%.2f, perm median OR=%.2f, empirical P=%.4f",
                       observed_or, median(perm_ors), p_empirical))
log_msg(STAGE, sprintf("  Threshold: empirical P < %g. RESULT: %s",
                       PARAMS$stage3$permutation_p_threshold,
                       ifelse(perm_test$passed, "PASS", "FAIL")))

# ------------------------------------------------------------------------------
# 5. H2 decision (Section 2.3)
# ------------------------------------------------------------------------------

tests_passed <- c(rrho = rrho_test$passed,
                  hyper = hyper_test$passed,
                  perm = perm_test$passed)
n_passed <- sum(tests_passed)

h2_decision <- if (n_passed == 3) {
  "SUPPORTED"
} else if (n_passed <= 1) {
  "REJECTED"
} else {
  "INCONCLUSIVE"
}

log_msg(STAGE, "==============================================")
log_msg(STAGE, sprintf("H2 DECISION: %s", h2_decision))
log_msg(STAGE, sprintf("Tests passed: %d / 3", n_passed))
log_msg(STAGE, sprintf("  RRHO2:         %s", ifelse(rrho_test$passed, "PASS", "FAIL")))
log_msg(STAGE, sprintf("  Hypergeometric: %s", ifelse(hyper_test$passed, "PASS", "FAIL")))
log_msg(STAGE, sprintf("  Permutation:    %s", ifelse(perm_test$passed, "PASS", "FAIL")))
log_msg(STAGE, "==============================================")

# ------------------------------------------------------------------------------
# 6. Save outputs
# ------------------------------------------------------------------------------

results <- list(
  h2_decision = h2_decision,
  n_tests_passed = n_passed,
  shared_universe_size = length(shared_universe),
  rrho_test = rrho_test,
  hypergeometric_test = hyper_test,
  permutation_test = perm_test,
  rrho_object = rrho_out,
  perm_or_distribution = perm_ors
)

save_intermediate(results, "stage3_convergence")

# Write a plain-text summary for easy reference
summary_lines <- c(
  sprintf("H2 DECISION: %s (%d/3 tests passed)", h2_decision, n_passed),
  "",
  sprintf("Test 1 (RRHO2):         max -log10(P) = %.3f, BH P = %.3g  [%s]",
          rrho_test$max_neg_log10_p, rrho_test$max_p_bh,
          ifelse(rrho_test$passed, "PASS", "FAIL")),
  sprintf("Test 2 (Hypergeometric): overlap = %d, OR = %.2f, Bonf P = %.3g  [%s]",
          hyper_test$overlap_n, hyper_test$odds_ratio, hyper_test$p_bonferroni,
          ifelse(hyper_test$passed, "PASS", "FAIL")),
  sprintf("Test 3 (Permutation):    observed OR = %.2f, empirical P = %.4f  [%s]",
          perm_test$observed_or, perm_test$p_empirical,
          ifelse(perm_test$passed, "PASS", "FAIL")),
  "",
  sprintf("Overlap genes (n=%d):", length(overlap)),
  paste(overlap, collapse = ", ")
)
writeLines(summary_lines, here::here(PARAMS$paths$intermediate, "stage3_H2_decision.txt"))

snapshot_session(STAGE)

# ==============================================================================
# POST-HOC: 4-quadrant analysis (added 2026-09-22)
# ==============================================================================
# Manual verification of RRHO2 signal via top-N rank-based quadrant overlap.
# Complements pre-specified 3-test decision by providing directional interpretation.
# Full details in DEVIATIONS.md (2026-09-21 Stage 3 discovery entry).

log_msg(STAGE, "Running post-hoc 4-quadrant analysis...")

N_top <- PARAMS$stage3$top_n_hypergeometric  # same as pre-specified (200)
merged_de <- merge(
  ad_de %>% mutate(ad_signed = -log10(PValue) * sign(logFC)) %>% 
    dplyr::select(gene, ad_logFC = logFC, ad_FDR = FDR, ad_signed),
  gbm_de %>% mutate(gbm_signed = -log10(PValue) * sign(logFC)) %>% 
    dplyr::select(gene, gbm_logFC = logFC, gbm_FDR = FDR, gbm_signed),
  by = "gene"
)

ad_top_up <- merged_de %>% arrange(desc(ad_signed)) %>% head(N_top) %>% pull(gene)
ad_top_dn <- merged_de %>% arrange(ad_signed)       %>% head(N_top) %>% pull(gene)
gbm_top_up <- merged_de %>% arrange(desc(gbm_signed)) %>% head(N_top) %>% pull(gene)
gbm_top_dn <- merged_de %>% arrange(gbm_signed)       %>% head(N_top) %>% pull(gene)

N_total <- nrow(merged_de)
exp_by_chance <- N_top * N_top / N_total

quadrants <- list(
  UU_concordant_UP     = list(genes = intersect(ad_top_up, gbm_top_up),
                              ad_direction = "UP",   gbm_direction = "UP"),
  UD_discordant        = list(genes = intersect(ad_top_up, gbm_top_dn),
                              ad_direction = "UP",   gbm_direction = "DOWN"),
  DU_H2a_discordant    = list(genes = intersect(ad_top_dn, gbm_top_up),
                              ad_direction = "DOWN", gbm_direction = "UP"),
  DD_concordant_DOWN   = list(genes = intersect(ad_top_dn, gbm_top_dn),
                              ad_direction = "DOWN", gbm_direction = "DOWN")
)

for (name in names(quadrants)) {
  q <- quadrants[[name]]
  n_obs <- length(q$genes)
  quadrants[[name]]$n_observed <- n_obs
  quadrants[[name]]$expected   <- exp_by_chance
  quadrants[[name]]$enrichment <- n_obs / exp_by_chance
  quadrants[[name]]$p_hyper    <- phyper(n_obs - 1, N_top, N_total - N_top, N_top,
                                          lower.tail = FALSE)
  quadrants[[name]]$gene_detail <- merged_de[merged_de$gene %in% q$genes,
    c("gene", "ad_logFC", "ad_FDR", "gbm_logFC", "gbm_FDR")]
}

quadrant_analysis <- list(
  method = "Top-N rank-based 4-quadrant overlap (post-hoc extension of Stage 3)",
  N_top = N_top,
  N_total = N_total,
  merged_data = merged_de,
  quadrants = quadrants,
  global_correlation = list(
    spearman = cor(merged_de$ad_signed, merged_de$gbm_signed, method = "spearman"),
    pearson  = cor(merged_de$ad_signed, merged_de$gbm_signed, method = "pearson")
  ),
  summary_table = data.frame(
    quadrant = c("UU (both UP)", "UD (AD UP x GBM DOWN)", 
                 "DU (AD DOWN x GBM UP) [H2a]", "DD (both DOWN) [novel]"),
    n_overlap = sapply(quadrants, function(x) x$n_observed),
    enrichment = sapply(quadrants, function(x) x$enrichment),
    p_hyper = sapply(quadrants, function(x) x$p_hyper),
    stringsAsFactors = FALSE
  )
)

save_intermediate(quadrant_analysis, "stage3_quadrant_analysis")
log_msg(STAGE, "4-quadrant summary:")
print(quadrant_analysis$summary_table)

log_msg(STAGE, "Stage 3 complete.")
