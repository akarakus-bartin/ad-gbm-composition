# ==============================================================================
# 04_stage1_ad_signature.R
# Purpose: Stage 1 - define AD-vulnerable RORB+ excitatory neuron signature
# from Leng et al. 2021 snRNA-seq data (GSE147528).
# Pre-specified in Analysis Plan Section 6.1.
# Outputs:
#   - ad_vulnerable_de_table (full edgeR results)
#   - ad_vulnerable_signature (ranked gene list, down-regulated in AD)
# ==============================================================================

source(here::here("R", "00_setup.R"))
STAGE <- "stage1"

leng <- load_intermediate("leng_seurat")

# ------------------------------------------------------------------------------
# 1. Subset: entorhinal cortex excitatory neurons only
# ------------------------------------------------------------------------------

# Region filter (SFG reserved for sensitivity per Section 8.1)
# NOTE: adapt column name to actual Leng metadata; typical column: "brain_region" or "region"
region_col <- if ("brain_region" %in% colnames(leng@meta.data)) "brain_region" else "region"
region_values <- leng[[region_col, drop = TRUE]]
ec_cells <- names(region_values)[grepl(PARAMS$stage1$brain_region, tolower(region_values))]
leng_ec <- subset(leng, cells = ec_cells)
log_msg(STAGE, sprintf("EC cells: %d", ncol(leng_ec)))

# Cell-type filter: excitatory neurons only
# NOTE: adapt column name to actual Leng metadata; typical column: "cell_type" or "celltype"
ct_col <- if ("cell_type" %in% colnames(leng_ec@meta.data)) "cell_type" else "celltype"
ct_values <- leng_ec[[ct_col, drop = TRUE]]
exc_cells <- names(ct_values)[grepl("excitatory|ExN|Exc", ct_values, ignore.case = TRUE)]
leng_exc <- subset(leng_ec, cells = exc_cells)
log_msg(STAGE, sprintf("EC excitatory neurons: %d", ncol(leng_exc)))

# ------------------------------------------------------------------------------
# 2. Define RORB-high subtype
# ------------------------------------------------------------------------------

if (!("RORB" %in% rownames(leng_exc))) {
  stop("RORB not detected in Leng expression matrix. Check gene naming convention (Ensembl vs Symbol).")
}

rorb_expr <- FetchData(leng_exc, vars = "RORB")[, "RORB"]
rorb_threshold <- quantile(rorb_expr[rorb_expr > 0], PARAMS$stage1$rorb_percentile)
log_msg(STAGE, sprintf("RORB expression 75th percentile among RORB+ cells: %.3f", rorb_threshold))

leng_exc$rorb_high <- rorb_expr >= rorb_threshold
log_msg(STAGE, sprintf("RORB-high cells: %d (%.1f%%)",
                       sum(leng_exc$rorb_high), 100 * mean(leng_exc$rorb_high)))

leng_rorb <- subset(leng_exc, subset = rorb_high == TRUE)

# ------------------------------------------------------------------------------
# 3. Donor-level pseudobulk (sum counts per donor)
# ------------------------------------------------------------------------------

# NOTE: adapt donor column name
donor_col <- if ("donor_id" %in% colnames(leng_rorb@meta.data)) "donor_id" else "sample_id"

# Extract raw counts
counts <- GetAssayData(leng_rorb, assay = "RNA", slot = "counts")
donor <- leng_rorb[[donor_col, drop = TRUE]]

# Pseudobulk: sum per donor
donors <- unique(donor)
pseudobulk <- sapply(donors, function(d) {
  cells <- names(donor)[donor == d]
  rowSums(counts[, cells, drop = FALSE])
})
colnames(pseudobulk) <- donors

log_msg(STAGE, sprintf("Pseudobulk matrix: %d genes x %d donors", nrow(pseudobulk), ncol(pseudobulk)))

# Donor metadata (one row per donor)
donor_meta <- leng_rorb@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  group_by(!!sym(donor_col)) %>%
  summarise(across(any_of(c("braak_stage", "diagnosis", "age", "sex")),
                   ~ first(.x)),
            n_cells = n(), .groups = "drop") %>%
  as.data.frame()
donor_meta <- donor_meta[match(colnames(pseudobulk), donor_meta[[donor_col]]), ]

# ------------------------------------------------------------------------------
# 4. Group assignment: Braak stage 0 (control) vs Braak stage VI (advanced AD)
# ------------------------------------------------------------------------------

braak_col <- "braak_stage"
donor_meta$group <- case_when(
  donor_meta[[braak_col]] %in% PARAMS$stage1$control_group_braak ~ "Control",
  donor_meta[[braak_col]] %in% PARAMS$stage1$ad_group_braak ~ "AD",
  TRUE ~ NA_character_
)

# Keep only Control and AD donors
keep <- !is.na(donor_meta$group)
donor_meta <- donor_meta[keep, ]
pseudobulk <- pseudobulk[, donor_meta[[donor_col]]]

n_ctrl <- sum(donor_meta$group == "Control")
n_ad <- sum(donor_meta$group == "AD")
log_msg(STAGE, sprintf("Donors: %d Control, %d AD", n_ctrl, n_ad))

assert_that(n_ctrl >= PARAMS$stage1$min_donors_per_group,
            sprintf("Too few Control donors (%d < %d)", n_ctrl, PARAMS$stage1$min_donors_per_group))
assert_that(n_ad >= PARAMS$stage1$min_donors_per_group,
            sprintf("Too few AD donors (%d < %d)", n_ad, PARAMS$stage1$min_donors_per_group))

# ------------------------------------------------------------------------------
# 5. edgeR differential expression
# ------------------------------------------------------------------------------

y <- DGEList(counts = pseudobulk, group = factor(donor_meta$group, levels = c("Control", "AD")))
keep_genes <- filterByExpr(y)
y <- y[keep_genes, , keep.lib.sizes = FALSE]
log_msg(STAGE, sprintf("Genes after filterByExpr: %d", nrow(y)))

y <- calcNormFactors(y, method = "TMM")

# Design: group + covariates (age, sex if available)
covariate_terms <- intersect(PARAMS$stage1$covariates, colnames(donor_meta))
if (length(covariate_terms) > 0) {
  # Only include covariates with variation
  covariate_terms <- covariate_terms[sapply(covariate_terms, function(v) length(unique(donor_meta[[v]])) > 1)]
}
design_formula <- as.formula(paste("~", paste(c("group", covariate_terms), collapse = " + ")))
design <- model.matrix(design_formula, data = donor_meta)
log_msg(STAGE, sprintf("Design formula: %s", deparse(design_formula)))

y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = "groupAD")

de_table <- topTags(qlf, n = Inf, sort.by = "PValue")$table %>%
  tibble::rownames_to_column("gene") %>%
  as_tibble()

log_msg(STAGE, sprintf("DE genes at FDR<%.2f: %d",
                       PARAMS$stage1$de_fdr_threshold,
                       sum(de_table$FDR < PARAMS$stage1$de_fdr_threshold)))

# ------------------------------------------------------------------------------
# 6. Define AD-vulnerable neuron down-regulated signature
# ------------------------------------------------------------------------------

signature <- de_table %>%
  filter(FDR < PARAMS$stage1$de_fdr_threshold,
         logFC < PARAMS$stage1$de_log2fc_threshold) %>%
  arrange(logFC) %>%
  mutate(rank = row_number())

log_msg(STAGE, sprintf("AD-vulnerable signature: %d down-regulated genes", nrow(signature)))
log_msg(STAGE, sprintf("Top 10: %s", paste(head(signature$gene, 10), collapse = ", ")))

# ------------------------------------------------------------------------------
# 7. Save outputs
# ------------------------------------------------------------------------------

save_intermediate(list(de_table = de_table,
                       signature = signature,
                       n_donors = c(control = n_ctrl, ad = n_ad),
                       rorb_threshold = rorb_threshold),
                  "stage1_ad_vulnerable")

# Also write the signature as a TSV for easy inspection
write_tsv(signature, here::here(PARAMS$paths$intermediate, "stage1_signature.tsv"))

snapshot_session(STAGE)
log_msg(STAGE, "Stage 1 complete.")
