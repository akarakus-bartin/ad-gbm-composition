# ==============================================================================
# 04_stage1_ad_signature.R
# Purpose: Stage 1 — define AD-vulnerable RORB+ excitatory neuron DE signature
# from Leng et al. 2021 snRNA-seq data (GSE147528) via pseudobulk edgeR.
#
# Pre-specified in Analysis Plan Section 6.1.
# Deviations documented in DEVIATIONS.md (2026-09-21 entry).
#
# Analysis design:
#   Primary:     Braak 0 vs 2 (RORB+ depletion peak per Leng Fig 2c)
#   Sensitivity: Braak 0 vs 6 (plan-conformant, late-AD control)
#
# Two models:
#   Model A (pooled):     ~ 0 + braak_group + subcluster
#                         Pooled DE across s1/s2/s4, subcluster as covariate.
#                         Primary output for H2 convergence testing.
#   Model B (per-cluster): DE within each subcluster separately.
#                          Consistency check across s1/s2/s4.
#
# Note: subcluster names contain ":" which is not R-syntactic. We sanitise
# via make.names() before building the design matrix so limma::makeContrasts
# accepts the column names.
# ==============================================================================

source(here::here("R", "00_setup.R"))
STAGE <- "stage1"

# ------------------------------------------------------------------------------
# Helper 1: aggregate SCE to pseudobulk (donor x subcluster) with min-cell filter
# ------------------------------------------------------------------------------
build_leng_pseudobulk <- function(sce, min_cells = 10) {
  group_id <- paste(as.character(sce$SampleID),
                    as.character(sce$subclusterAssignment), sep = "__")
  agg <- scuttle::aggregateAcrossCells(sce, ids = group_id,
                                        statistics = "sum",
                                        use.assay.type = "counts")
  n_cells <- as.integer(table(group_id)[colnames(agg)])
  agg$n_cells <- n_cells
  keep <- agg$n_cells >= min_cells
  log_msg(STAGE, sprintf("Pseudobulk samples: %d total, %d pass min_cells=%d",
                         ncol(agg), sum(keep), min_cells))
  agg <- agg[, keep]
  sample_meta <- data.frame(
    sample_id           = colnames(agg),
    donor               = as.character(agg$SampleID),
    subcluster_original = as.character(agg$subclusterAssignment),
    subcluster          = make.names(as.character(agg$subclusterAssignment)),
    braak_stage         = as.character(agg$BraakStage),
    n_cells             = agg$n_cells,
    stringsAsFactors    = FALSE
  )
  list(counts = SummarizedExperiment::assay(agg, "counts"), meta = sample_meta)
}

# ------------------------------------------------------------------------------
# Helper 2: Model A pooled DE (Braak_test vs Braak_ref, subcluster covariate)
# ------------------------------------------------------------------------------
run_edger_pooled <- function(counts, meta, braak_ref, braak_test) {
  keep <- meta$braak_stage %in% c(braak_ref, braak_test)
  counts_sub <- counts[, keep]
  meta_sub <- meta[keep, , drop = FALSE]
  meta_sub$braak_group <- factor(paste0("Braak", meta_sub$braak_stage),
                                  levels = c(paste0("Braak", braak_ref),
                                             paste0("Braak", braak_test)))
  meta_sub$subcluster <- factor(meta_sub$subcluster)
  design <- model.matrix(~ 0 + braak_group + subcluster, data = meta_sub)
  y <- edgeR::DGEList(counts = counts_sub, samples = meta_sub)
  keep_g <- edgeR::filterByExpr(y, design = design)
  y <- y[keep_g, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y, method = "TMM")
  y <- edgeR::estimateDisp(y, design = design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design = design, robust = TRUE)
  contrast_str <- sprintf("braak_groupBraak%s - braak_groupBraak%s",
                          braak_test, braak_ref)
  contrast <- limma::makeContrasts(contrasts = contrast_str, levels = design)
  qlf <- edgeR::glmQLFTest(fit, contrast = contrast)
  tt <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
  tt$gene <- rownames(tt)
  tt <- tt[, c("gene", "logFC", "logCPM", "F", "PValue", "FDR")]
  tt <- tt[order(tt$PValue), ]
  log_msg(STAGE, sprintf("Model A [Braak %s vs %s]: %d genes, %d FDR<0.05",
                         braak_ref, braak_test, nrow(tt),
                         sum(tt$FDR < 0.05, na.rm = TRUE)))
  list(de_table = tt, n_samples = ncol(counts_sub),
       n_donors_ref  = length(unique(meta_sub$donor[meta_sub$braak_stage == braak_ref])),
       n_donors_test = length(unique(meta_sub$donor[meta_sub$braak_stage == braak_test])),
       contrast = contrast_str)
}

# ------------------------------------------------------------------------------
# Helper 3: Model B per-subcluster DE
# ------------------------------------------------------------------------------
run_edger_per_subcluster <- function(counts, meta, braak_ref, braak_test, subclusters) {
  results <- list()
  for (sc in subclusters) {
    keep <- meta$subcluster == sc & meta$braak_stage %in% c(braak_ref, braak_test)
    if (sum(keep) < 4) {
      log_msg(STAGE, sprintf("  %s: skipped (only %d samples)", sc, sum(keep)))
      next
    }
    counts_sub <- counts[, keep]
    meta_sub <- meta[keep, , drop = FALSE]
    meta_sub$braak_group <- factor(paste0("Braak", meta_sub$braak_stage),
                                    levels = c(paste0("Braak", braak_ref),
                                               paste0("Braak", braak_test)))
    design <- model.matrix(~ braak_group, data = meta_sub)
    y <- edgeR::DGEList(counts = counts_sub, samples = meta_sub)
    keep_g <- edgeR::filterByExpr(y, design = design)
    y <- y[keep_g, , keep.lib.sizes = FALSE]
    y <- edgeR::calcNormFactors(y, method = "TMM")
    y <- edgeR::estimateDisp(y, design = design, robust = TRUE)
    fit <- edgeR::glmQLFit(y, design = design, robust = TRUE)
    qlf <- edgeR::glmQLFTest(fit, coef = 2)
    tt <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
    tt$gene <- rownames(tt)
    tt <- tt[, c("gene", "logFC", "logCPM", "F", "PValue", "FDR")]
    tt <- tt[order(tt$PValue), ]
    log_msg(STAGE, sprintf("  %s (n=%d): %d genes, %d FDR<0.05",
                            sc, ncol(counts_sub), nrow(tt),
                            sum(tt$FDR < 0.05, na.rm = TRUE)))
    results[[sc]] <- tt
  }
  results
}

# ==============================================================================
# Main pipeline
# ==============================================================================

leng <- load_intermediate("leng_processed")
assert_that(!is.null(leng$sce_vulnerable), "leng_processed missing sce_vulnerable")
log_msg(STAGE, sprintf("Loaded Leng vulnerable subset: %d cells x %d genes",
                       ncol(leng$sce_vulnerable), nrow(leng$sce_vulnerable)))

min_cells <- if (!is.null(PARAMS$leng$min_cells_per_pseudobulk)) PARAMS$leng$min_cells_per_pseudobulk else 10
pb <- build_leng_pseudobulk(leng$sce_vulnerable, min_cells = min_cells)
log_msg(STAGE, sprintf("Pseudobulk matrix: %d genes x %d samples",
                       nrow(pb$counts), ncol(pb$counts)))
log_msg(STAGE, sprintf("Braak distribution: %s",
                       paste(names(table(pb$meta$braak_stage)),
                             table(pb$meta$braak_stage), sep = "=", collapse = ", ")))
save_intermediate(pb, "leng_pseudobulk")

log_msg(STAGE, "=== Model A (pooled, subcluster covariate) ===")
log_msg(STAGE, "Primary contrast: Braak 0 vs 2")
deg_primary <- run_edger_pooled(pb$counts, pb$meta, "0", "2")
save_intermediate(deg_primary, "leng_deg_primary")

log_msg(STAGE, "Sensitivity contrast: Braak 0 vs 6")
deg_sensitivity <- run_edger_pooled(pb$counts, pb$meta, "0", "6")
save_intermediate(deg_sensitivity, "leng_deg_sensitivity")

log_msg(STAGE, "=== Model B (per-subcluster, Braak 0 vs 2) ===")
vuln_sc <- make.names(leng$vulnerable_subclusters)
deg_per_sc <- run_edger_per_subcluster(pb$counts, pb$meta, "0", "2", vuln_sc)
save_intermediate(deg_per_sc, "leng_deg_per_subcluster")

log_msg(STAGE, "Stage 1 complete.")
