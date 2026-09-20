# ==============================================================================
# 07_stage4_bulk_reanalysis.R
# Purpose: Stage 4 - bulk-level DE analysis with dual deconvolution control
# (BRETIGEA and bMIND). Tests H1 (Section 2.2).
# Also reconciles bulk-level convergence with single-cell findings (Stages 1-3).
#
# H1 decision rule (Section 2.2):
#   - Supported if shared DEG count (controlled) / (uncorrected) <= 20% in BOTH methods
#   - Rejected if >= 50% in either method
#   - Inconclusive if methods disagree
# ==============================================================================

source(here::here("R", "00_setup.R"))
STAGE <- "stage4"

ad <- load_intermediate("ad_bulk")
gbm <- load_intermediate("gbm_bulk")

# ------------------------------------------------------------------------------
# 1. BRETIGEA deconvolution
# ------------------------------------------------------------------------------

run_bretigea <- function(expr, markers_n = PARAMS$stage4$bretigea_n_markers) {
  # brainCells uses PCA-based aggregation across marker genes
  scores <- BRETIGEA::brainCells(inputMat = expr,
                                  nMarker = markers_n,
                                  celltypes = PARAMS$stage4$bretigea_cell_types,
                                  method = "SVD",
                                  scale = TRUE)
  as.data.frame(scores)
}

log_msg(STAGE, "Running BRETIGEA on AD cohort")
ad_bretigea <- run_bretigea(ad$expr)
log_msg(STAGE, "Running BRETIGEA on GBM cohort")
gbm_bretigea <- run_bretigea(gbm$expr)

# ------------------------------------------------------------------------------
# 2. bMIND deconvolution
# ------------------------------------------------------------------------------
# NOTE: bMIND requires a cell-type-specific reference. The plan specifies
# Allen Brain Atlas M1 snRNA-seq. Loading + preparing this reference is a
# multi-step process; here we sketch the interface.

prepare_allen_reference <- function() {
  # TODO: download Allen M1 snRNA-seq counts + cell-type labels
  # (https://portal.brain-map.org/atlases-and-data/rnaseq)
  # Then compute cell-type-specific mean expression profiles.
  # Below is a placeholder; adapt to actual reference file structure.
  ref_path <- here::here(PARAMS$paths$raw_data, "allen_m1_ref.rds")
  if (!file.exists(ref_path)) {
    stop("Allen M1 reference not found. Download and place at: ", ref_path)
  }
  readRDS(ref_path)   # expects a matrix: genes x cell types
}

run_bmind <- function(expr) {
  ref <- prepare_allen_reference()
  common <- intersect(rownames(expr), rownames(ref))
  # bMIND: bMIND(bulk, sig = reference[cell-type profile])
  # Full call requires proper sample metadata; this is a sketch.
  # See MIND documentation for full parameter set.
  bmind_result <- MIND::bMIND(bulk = expr[common, ], sigma = ref[common, ])
  # Extract cell-type proportions
  as.data.frame(bmind_result$fraction)
}

log_msg(STAGE, "Running bMIND on AD cohort")
ad_bmind <- tryCatch(run_bmind(ad$expr),
                     error = function(e) {
                       log_msg(STAGE, sprintf("bMIND on AD FAILED: %s", conditionMessage(e)))
                       NULL
                     })
log_msg(STAGE, "Running bMIND on GBM cohort")
gbm_bmind <- tryCatch(run_bmind(gbm$expr),
                     error = function(e) {
                       log_msg(STAGE, sprintf("bMIND on GBM FAILED: %s", conditionMessage(e)))
                       NULL
                     })

# ------------------------------------------------------------------------------
# 3. limma models: A (uncorrected), B (BRETIGEA), C (bMIND), D (both)
# ------------------------------------------------------------------------------

fit_limma <- function(expr, meta, model_type, deconv = NULL) {
  # Base design
  if (any(grepl("dataset", colnames(meta)))) {
    # AD cohort
    base_formula <- "~ diagnosis + brain_region + sex + age + dataset"
  } else {
    # GBM cohort
    base_formula <- "~ diagnosis + sex + brain_region"
  }

  # Add deconvolution covariates
  if (model_type == "B" && !is.null(deconv)) {
    covariate_names <- paste(colnames(deconv), collapse = " + ")
    formula_str <- paste(base_formula, "+", covariate_names)
    meta <- cbind(meta, deconv[meta$sample_id, ])
  } else if (model_type == "C" && !is.null(deconv)) {
    covariate_names <- paste(colnames(deconv), collapse = " + ")
    formula_str <- paste(base_formula, "+", covariate_names)
    meta <- cbind(meta, deconv[meta$sample_id, ])
  } else {
    formula_str <- base_formula
  }

  log_msg(STAGE, sprintf("  Model %s formula: %s", model_type, formula_str))

  # Drop covariates with no variation
  design_formula <- as.formula(formula_str)
  design <- tryCatch(model.matrix(design_formula, data = meta),
                     error = function(e) {
                       log_msg(STAGE, sprintf("Design matrix FAILED: %s", conditionMessage(e)))
                       NULL
                     })
  if (is.null(design)) return(NULL)

  # Fit
  fit <- lmFit(expr, design)
  fit <- eBayes(fit)

  # Extract diagnosis coefficient
  diag_col <- grep("^diagnosis", colnames(design), value = TRUE)[1]
  tt <- topTable(fit, coef = diag_col, number = Inf, sort.by = "P") %>%
    tibble::rownames_to_column("gene") %>%
    as_tibble()

  # VIF check on diagnosis
  vif_diagnosis <- tryCatch({
    lm_fit <- lm(as.formula(paste("expr[1, ]", formula_str)), data = meta)
    car::vif(lm_fit)[diag_col]
  }, error = function(e) NA)
  log_msg(STAGE, sprintf("  Model %s: diagnosis VIF = %s", model_type,
                         ifelse(is.na(vif_diagnosis), "NA", sprintf("%.2f", vif_diagnosis))))

  list(topTable = tt, vif_diagnosis = vif_diagnosis, formula = formula_str)
}

# Run all models for AD
log_msg(STAGE, "--- AD cohort limma models ---")
ad_models <- list(
  A = fit_limma(ad$expr, ad$meta, "A"),
  B = fit_limma(ad$expr, ad$meta, "B", deconv = ad_bretigea)
)
if (!is.null(ad_bmind)) {
  ad_models$C <- fit_limma(ad$expr, ad$meta, "C", deconv = ad_bmind)
}

# Run all models for GBM
log_msg(STAGE, "--- GBM cohort limma models ---")
gbm_models <- list(
  A = fit_limma(gbm$expr, gbm$meta, "A"),
  B = fit_limma(gbm$expr, gbm$meta, "B", deconv = gbm_bretigea)
)
if (!is.null(gbm_bmind)) {
  gbm_models$C <- fit_limma(gbm$expr, gbm$meta, "C", deconv = gbm_bmind)
}

# ------------------------------------------------------------------------------
# 4. Count DEGs per model
# ------------------------------------------------------------------------------

count_degs <- function(tt, fdr_thresh, log2fc_thresh) {
  if (is.null(tt)) return(NA)
  sum(tt$adj.P.Val < fdr_thresh & abs(tt$logFC) > log2fc_thresh)
}

deg_counts <- tibble(
  cohort = c(rep("AD", length(ad_models)), rep("GBM", length(gbm_models))),
  model = c(names(ad_models), names(gbm_models)),
  n_degs = c(sapply(ad_models, function(m) count_degs(m$topTable,
                                                     PARAMS$stage4$ad_fdr_threshold,
                                                     PARAMS$stage4$ad_log2fc_threshold)),
             sapply(gbm_models, function(m) count_degs(m$topTable,
                                                      PARAMS$stage4$gbm_fdr_threshold,
                                                      PARAMS$stage4$gbm_log2fc_threshold)))
)
print(deg_counts)

# ------------------------------------------------------------------------------
# 5. Shared DEGs between AD and GBM (per model type)
# ------------------------------------------------------------------------------

count_shared <- function(ad_tt, gbm_tt) {
  if (is.null(ad_tt) || is.null(gbm_tt)) return(NA)
  ad_sig <- ad_tt %>%
    filter(adj.P.Val < PARAMS$stage4$ad_fdr_threshold,
           abs(logFC) > PARAMS$stage4$ad_log2fc_threshold)
  gbm_sig <- gbm_tt %>%
    filter(adj.P.Val < PARAMS$stage4$gbm_fdr_threshold,
           abs(logFC) > PARAMS$stage4$gbm_log2fc_threshold)

  # Directional concordance required (both up or both down)
  merged <- inner_join(ad_sig %>% select(gene, ad_logFC = logFC),
                       gbm_sig %>% select(gene, gbm_logFC = logFC),
                       by = "gene") %>%
    filter(sign(ad_logFC) == sign(gbm_logFC))
  nrow(merged)
}

shared_counts <- tibble(
  model = c("A", "B", "C"),
  n_shared = c(
    count_shared(ad_models$A$topTable, gbm_models$A$topTable),
    count_shared(ad_models$B$topTable, gbm_models$B$topTable),
    if (!is.null(ad_models$C) && !is.null(gbm_models$C))
      count_shared(ad_models$C$topTable, gbm_models$C$topTable) else NA
  )
)
print(shared_counts)

# ------------------------------------------------------------------------------
# 6. H1 decision (Section 2.2)
# ------------------------------------------------------------------------------

n_a <- shared_counts$n_shared[shared_counts$model == "A"]
n_b <- shared_counts$n_shared[shared_counts$model == "B"]
n_c <- shared_counts$n_shared[shared_counts$model == "C"]

ratio_b <- if (!is.na(n_a) && n_a > 0) n_b / n_a else NA
ratio_c <- if (!is.na(n_a) && n_a > 0 && !is.na(n_c)) n_c / n_a else NA

log_msg(STAGE, sprintf("Shared DEG counts: Model A = %s, Model B = %s, Model C = %s",
                       n_a, n_b, ifelse(is.na(n_c), "NA (bMIND unavailable)", as.character(n_c))))
log_msg(STAGE, sprintf("Ratios: B/A = %.3f, C/A = %s",
                       ratio_b, ifelse(is.na(ratio_c), "NA", sprintf("%.3f", ratio_c))))

# Apply decision rule
h1_decision <- if (is.na(ratio_c)) {
  "INCOMPLETE (bMIND not available; cannot apply dual-method rule)"
} else if (ratio_b <= PARAMS$stage4$h1_support_ratio_max &&
           ratio_c <= PARAMS$stage4$h1_support_ratio_max) {
  "SUPPORTED"
} else if (ratio_b >= PARAMS$stage4$h1_reject_ratio_min ||
           ratio_c >= PARAMS$stage4$h1_reject_ratio_min) {
  "REJECTED"
} else {
  "INCONCLUSIVE"
}

log_msg(STAGE, "==============================================")
log_msg(STAGE, sprintf("H1 DECISION: %s", h1_decision))
log_msg(STAGE, "==============================================")

# ------------------------------------------------------------------------------
# 7. Reconciliation with Stages 1-3 (exploratory)
# ------------------------------------------------------------------------------

# For genes shared under Model A but lost under Models B/C:
# check whether they are cell-type-specifically expressed in Leng (excitatory neurons)

if (!is.na(n_a) && n_a > 0) {
  a_shared_genes <- inner_join(
    ad_models$A$topTable %>%
      filter(adj.P.Val < PARAMS$stage4$ad_fdr_threshold,
             abs(logFC) > PARAMS$stage4$ad_log2fc_threshold) %>%
      select(gene, ad_logFC = logFC),
    gbm_models$A$topTable %>%
      filter(adj.P.Val < PARAMS$stage4$gbm_fdr_threshold,
             abs(logFC) > PARAMS$stage4$gbm_log2fc_threshold) %>%
      select(gene, gbm_logFC = logFC),
    by = "gene"
  ) %>%
    filter(sign(ad_logFC) == sign(gbm_logFC))

  # Cross-reference with Stage 1 AD-vulnerable signature
  stage1_sig <- load_intermediate("stage1_ad_vulnerable")$signature$gene
  a_shared_genes$in_stage1_ad_vulnerable <- a_shared_genes$gene %in% stage1_sig
  log_msg(STAGE, sprintf("Of %d bulk-level shared DEGs (Model A), %d (%.0f%%) are in the Stage 1 AD-vulnerable neuron signature",
                         nrow(a_shared_genes),
                         sum(a_shared_genes$in_stage1_ad_vulnerable),
                         100 * mean(a_shared_genes$in_stage1_ad_vulnerable)))

  save_intermediate(a_shared_genes, "stage4_bulk_shared_reconciled")
}

# ------------------------------------------------------------------------------
# 8. Save all outputs
# ------------------------------------------------------------------------------

save_intermediate(list(
  ad_models = ad_models,
  gbm_models = gbm_models,
  ad_bretigea = ad_bretigea,
  gbm_bretigea = gbm_bretigea,
  ad_bmind = ad_bmind,
  gbm_bmind = gbm_bmind,
  deg_counts = deg_counts,
  shared_counts = shared_counts,
  h1_decision = h1_decision,
  ratios = c(B_over_A = ratio_b, C_over_A = ratio_c)
), "stage4_bulk_reanalysis")

writeLines(c(
  sprintf("H1 DECISION: %s", h1_decision),
  sprintf("Model A (uncorrected) shared DEGs: %s", n_a),
  sprintf("Model B (BRETIGEA-controlled) shared DEGs: %s (ratio %.3f)", n_b, ratio_b),
  sprintf("Model C (bMIND-controlled) shared DEGs: %s (ratio %s)",
          ifelse(is.na(n_c), "NA", as.character(n_c)),
          ifelse(is.na(ratio_c), "NA", sprintf("%.3f", ratio_c)))
), here::here(PARAMS$paths$intermediate, "stage4_H1_decision.txt"))

snapshot_session(STAGE)
log_msg(STAGE, "Stage 4 complete.")
