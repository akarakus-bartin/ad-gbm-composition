# ==============================================================================
# 08_sensitivity.R
# Purpose: run all sensitivity analyses specified in Analysis Plan Section 8.
# Each section maps directly to Section 8.1-8.4 of the plan.
# Results are labelled and saved separately so they cannot be confused with
# the primary analysis.
# ==============================================================================

source(here::here("R", "00_setup.R"))
STAGE <- "sensitivity"

leng <- load_intermediate("leng_seurat")
neftel <- load_intermediate("neftel_seurat")
ad_bulk <- load_intermediate("ad_bulk")
gbm_bulk <- load_intermediate("gbm_bulk")
stage1 <- load_intermediate("stage1_ad_vulnerable")
stage2 <- load_intermediate("stage2_gbm_neural_mimicry")
stage3 <- load_intermediate("stage3_convergence")
stage4 <- load_intermediate("stage4_bulk_reanalysis")

sensitivity_results <- list()

# ==============================================================================
# 8.1 Sensitivity of the AD-vulnerable neuron signature
# ==============================================================================

log_msg(STAGE, "8.1.a Intermediate Braak stage II as AD group")
# Re-run Stage 1 with Braak II instead of Braak VI
# Load the Stage 1 script logic and swap PARAMS$stage1$ad_group_braak = c(2)
# Simplest: source Stage 1 in a modified environment
# NOTE: implement as a function refactor in a real project; here we sketch.

sensitivity_results$s8_1a_braak2 <- list(
  description = "Stage 1 repeated with Braak stage II as AD group",
  status = "TODO: refactor Stage 1 as a function and re-invoke with Braak = c(2)",
  planned_output = "AD-vulnerable signature (Braak II) overlap with primary (Braak VI) signature"
)

log_msg(STAGE, "8.1.b RORB threshold at 66th and 90th percentiles")
sensitivity_results$s8_1b_rorb <- list(
  description = "Stage 1 repeated with RORB threshold at 66th and 90th percentiles",
  status = "TODO: re-run with PARAMS$stage1$rorb_percentile in c(0.66, 0.90)",
  planned_output = "Overlap of signatures across three thresholds; concordance measure"
)

log_msg(STAGE, "8.1.c Include SFG as second brain region")
sensitivity_results$s8_1c_sfg <- list(
  description = "Stage 1 repeated including superior frontal gyrus samples",
  status = "TODO: re-run with region filter widened to include SFG",
  planned_output = "SFG-derived signature vs EC-derived signature comparison"
)

# ==============================================================================
# 8.2 Sensitivity of the GBM neural-mimicry signature
# ==============================================================================

log_msg(STAGE, "8.2.a NPC-like only (excluding OPC-like)")
sensitivity_results$s8_2a_npc_only <- list(
  description = "Stage 2 repeated with only NPC-like as neural-lineage group",
  status = "TODO: re-run with PARAMS$stage2$neural_lineage_states = c('NPClike')",
  planned_output = "Signature-signature overlap; Test 3 rerun if applicable"
)

log_msg(STAGE, "8.2.b Suva 2014 developmental hierarchy classifier")
sensitivity_results$s8_2b_suva <- list(
  description = "Stage 2 repeated using Suva et al. 2014 developmental hierarchy states",
  status = "TODO: implement Suva score computation; re-classify cells; re-run DE",
  planned_output = "Alternative state assignment; effect on neural-mimicry signature"
)

# ==============================================================================
# 8.3 Sensitivity of the convergence test
# ==============================================================================

# 8.3.a Test 2 with different top-N thresholds
log_msg(STAGE, "8.3.a Test 2 with top-100 and top-500 genes")

test_hypergeometric <- function(N, ad_de, gbm_de, universe) {
  ad_top <- ad_de %>% arrange(logFC) %>% slice_head(n = N) %>% pull(gene)
  gbm_top <- gbm_de %>% arrange(desc(logFC)) %>% slice_head(n = N) %>% pull(gene)
  overlap <- length(intersect(ad_top, gbm_top))
  p_hyper <- phyper(q = overlap - 1, m = N, n = length(universe) - N, k = N, lower.tail = FALSE)
  a <- overlap; b <- N - a; c <- N - a; d <- length(universe) - a - b - c
  or <- (a * d) / max(1, b * c)
  list(N = N, overlap = overlap, odds_ratio = or, p_raw = p_hyper)
}

shared_universe <- intersect(stage1$de_table$gene, stage2$de_table$gene)
ad_de_shared <- stage1$de_table %>% filter(gene %in% shared_universe)
gbm_de_shared <- stage2$de_table %>% filter(gene %in% shared_universe)

sensitivity_results$s8_3a <- list(
  n100 = test_hypergeometric(100, ad_de_shared, gbm_de_shared, shared_universe),
  n200_primary = test_hypergeometric(200, ad_de_shared, gbm_de_shared, shared_universe),
  n500 = test_hypergeometric(500, ad_de_shared, gbm_de_shared, shared_universe)
)
log_msg(STAGE, sprintf("Test 2 by N: 100 -> OR=%.2f (P=%.3g); 200 -> OR=%.2f (P=%.3g); 500 -> OR=%.2f (P=%.3g)",
                       sensitivity_results$s8_3a$n100$odds_ratio, sensitivity_results$s8_3a$n100$p_raw,
                       sensitivity_results$s8_3a$n200_primary$odds_ratio, sensitivity_results$s8_3a$n200_primary$p_raw,
                       sensitivity_results$s8_3a$n500$odds_ratio, sensitivity_results$s8_3a$n500$p_raw))

# 8.3.b Test 3 with 10,000 permutations
log_msg(STAGE, "8.3.b Test 3 with 10,000 permutations (this may take a while)")
N <- PARAMS$stage3$top_n_hypergeometric
set.seed(PARAMS$seed)
perm_ors_10k <- replicate(10000, {
  set1 <- sample(shared_universe, N); set2 <- sample(shared_universe, N)
  a <- length(intersect(set1, set2)); b <- N - a; c <- N - a
  d <- length(shared_universe) - a - b - c
  (a * d) / max(1, b * c)
})
observed_or <- stage3$permutation_test$observed_or
p_empirical_10k <- (sum(perm_ors_10k >= observed_or) + 1) / (10001)
sensitivity_results$s8_3b_10k_permutations <- list(
  observed_or = observed_or,
  p_empirical_10k = p_empirical_10k,
  p_empirical_1k_primary = stage3$permutation_test$p_empirical
)
log_msg(STAGE, sprintf("Test 3 (10k perm): P=%.5f (vs 1k perm P=%.5f)",
                       p_empirical_10k, stage3$permutation_test$p_empirical))

# 8.3.c Restrict to protein-coding only
log_msg(STAGE, "8.3.c Restricting to protein-coding genes")
# NOTE: requires an annotation table (biomaRt or org.Hs.eg.db)
# Placeholder:
sensitivity_results$s8_3c_protein_coding <- list(
  description = "Repeat all 3 tests restricted to protein-coding genes",
  status = "TODO: filter shared_universe by gene_biotype == 'protein_coding' via biomaRt"
)

# ==============================================================================
# 8.4 Sensitivity of the bulk deconvolution result
# ==============================================================================

log_msg(STAGE, "8.4.a CIBERSORTx triangulation")
sensitivity_results$s8_4a_cibersortx <- list(
  description = "Repeat Stage 4 with CIBERSORTx as a third deconvolution method",
  status = "TODO: CIBERSORTx runs on a web service or via CLI with a signature matrix. Adapt to your setup."
)

log_msg(STAGE, "8.4.b Unweighted average of BRETIGEA and bMIND")
sensitivity_results$s8_4b_average <- list(
  description = "Repeat Stage 4 with averaged BRETIGEA + bMIND cell-type scores",
  status = "TODO: rescale both to z-scores per cell type, average, use as covariates in limma"
)

log_msg(STAGE, "8.4.c Exclude low-purity GBM samples")
if (!is.null(gbm_bulk$meta$tumor_purity)) {
  purity_ok <- gbm_bulk$meta$tumor_purity >= 0.6 | is.na(gbm_bulk$meta$tumor_purity)
  log_msg(STAGE, sprintf("Samples with purity >= 60%%: %d / %d",
                         sum(purity_ok, na.rm = TRUE), length(purity_ok)))
  sensitivity_results$s8_4c_purity <- list(
    description = "Stage 4 GBM analysis excluding samples with tumor purity < 60%",
    status = "TODO: re-run Stage 4 on subsetted GBM cohort",
    n_samples_retained = sum(purity_ok, na.rm = TRUE)
  )
} else {
  sensitivity_results$s8_4c_purity <- list(
    description = "Tumor purity data not available in metadata",
    status = "NEED: fetch tumor purity from ABSOLUTE (TCGA) or ESTIMATE algorithm"
  )
}

# ------------------------------------------------------------------------------
# Save sensitivity results
# ------------------------------------------------------------------------------

save_intermediate(sensitivity_results, "sensitivity_results")

snapshot_session(STAGE)
log_msg(STAGE, "Sensitivity analyses complete (with TODO markers where noted).")
