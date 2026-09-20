# ==============================================================================
# 05_stage2_gbm_signature.R
# Purpose: Stage 2 - characterise GBM neural-mimicry signature from Neftel et al.
# 2019 scRNA-seq data (GSE131928).
# Pre-specified in Analysis Plan Section 6.2.
# Contrast: (NPClike + OPClike) vs (AClike + MESlike)
# Outputs:
#   - gbm_neural_mimicry_de_table
#   - gbm_neural_mimicry_signature (ranked gene list, up-regulated in neural-lineage states)
# ==============================================================================

source(here::here("R", "00_setup.R"))
STAGE <- "stage2"

neftel <- load_intermediate("neftel_seurat")

# ------------------------------------------------------------------------------
# 1. Verify state assignment
# ------------------------------------------------------------------------------

state_col <- if ("cellular_tumor_state" %in% colnames(neftel@meta.data)) {
  "cellular_tumor_state"
} else if ("state" %in% colnames(neftel@meta.data)) {
  "state"
} else {
  stop("Neftel cellular state column not found. Adapt state_col in this script.")
}

state <- neftel[[state_col, drop = TRUE]]
log_msg(STAGE, "State distribution:")
print(table(state))

# Keep only cells in the 4 defined states (drop intermediates / unassigned)
keep_states <- c(PARAMS$stage2$neural_lineage_states, PARAMS$stage2$other_states)
keep_cells <- names(state)[state %in% keep_states]
neftel_f <- subset(neftel, cells = keep_cells)

# Contrast group
group <- ifelse(state[keep_cells] %in% PARAMS$stage2$neural_lineage_states,
                "NeuralLineage", "Other")
neftel_f$state_contrast <- group
log_msg(STAGE, sprintf("Neural-lineage (NPClike + OPClike): %d cells", sum(group == "NeuralLineage")))
log_msg(STAGE, sprintf("Other (AClike + MESlike): %d cells", sum(group == "Other")))

# ------------------------------------------------------------------------------
# 2. Pseudobulk by tumor x state_contrast
# ------------------------------------------------------------------------------

tumor_col <- if ("tumor_id" %in% colnames(neftel_f@meta.data)) "tumor_id" else "sample_id"
tumor <- neftel_f[[tumor_col, drop = TRUE]]

counts <- GetAssayData(neftel_f, assay = "RNA", slot = "counts")

# Build combinations of tumor x contrast
meta <- data.frame(cell = colnames(counts),
                   tumor = tumor,
                   contrast = neftel_f$state_contrast,
                   stringsAsFactors = FALSE)
meta$pb_id <- paste(meta$tumor, meta$contrast, sep = "_")

# Sum counts per pseudobulk group
pb_ids <- unique(meta$pb_id)
pseudobulk <- sapply(pb_ids, function(id) {
  cells <- meta$cell[meta$pb_id == id]
  rowSums(counts[, cells, drop = FALSE])
})
colnames(pseudobulk) <- pb_ids

# Pseudobulk-level metadata
pb_meta <- meta %>%
  distinct(pb_id, tumor, contrast) %>%
  mutate(n_cells = sapply(pb_id, function(id) sum(meta$pb_id == id))) %>%
  filter(n_cells >= PARAMS$stage2$min_cells_per_tumor_state)

# Drop pseudobulks with too few cells
pseudobulk <- pseudobulk[, pb_meta$pb_id]
log_msg(STAGE, sprintf("Pseudobulk matrix: %d genes x %d tumor-state groups", nrow(pseudobulk), ncol(pseudobulk)))
log_msg(STAGE, sprintf("Tumors with both contrast groups: %d",
                       length(intersect(pb_meta$tumor[pb_meta$contrast == "NeuralLineage"],
                                        pb_meta$tumor[pb_meta$contrast == "Other"]))))

# ------------------------------------------------------------------------------
# 3. edgeR DE with tumor as blocking factor
# ------------------------------------------------------------------------------

y <- DGEList(counts = pseudobulk)
keep_genes <- filterByExpr(y, group = pb_meta$contrast)
y <- y[keep_genes, , keep.lib.sizes = FALSE]
log_msg(STAGE, sprintf("Genes after filterByExpr: %d", nrow(y)))
y <- calcNormFactors(y, method = "TMM")

# Design: contrast + tumor (blocking)
design <- model.matrix(~ 0 + contrast + tumor, data = pb_meta)
colnames(design) <- make.names(colnames(design))

y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

# Contrast: NeuralLineage vs Other
contrast_matrix <- makeContrasts(
  NL_vs_Other = contrastNeuralLineage - contrastOther,
  levels = design
)
qlf <- glmQLFTest(fit, contrast = contrast_matrix)

de_table <- topTags(qlf, n = Inf, sort.by = "PValue")$table %>%
  tibble::rownames_to_column("gene") %>%
  as_tibble()

log_msg(STAGE, sprintf("DE genes at FDR<%.2f: %d",
                       PARAMS$stage2$de_fdr_threshold,
                       sum(de_table$FDR < PARAMS$stage2$de_fdr_threshold)))

# ------------------------------------------------------------------------------
# 4. Define GBM neural-mimicry up-regulated signature
# ------------------------------------------------------------------------------

signature <- de_table %>%
  filter(FDR < PARAMS$stage2$de_fdr_threshold,
         logFC > PARAMS$stage2$de_log2fc_threshold) %>%
  arrange(desc(logFC)) %>%
  mutate(rank = row_number())

log_msg(STAGE, sprintf("GBM neural-mimicry signature: %d up-regulated genes", nrow(signature)))
log_msg(STAGE, sprintf("Top 10: %s", paste(head(signature$gene, 10), collapse = ", ")))

# ------------------------------------------------------------------------------
# 5. Internal positive control: check presence of known presynaptic markers
# ------------------------------------------------------------------------------

presynaptic_markers <- c("SYN1", "SYT1", "STXBP1", "NAPB", "VAMP2", "SNAP25", "RAB3A", "SLC17A7")
present <- intersect(presynaptic_markers, signature$gene)
log_msg(STAGE, sprintf("Presynaptic markers found in signature (%d/%d): %s",
                       length(present), length(presynaptic_markers),
                       paste(present, collapse = ", ")))

# ------------------------------------------------------------------------------
# 6. Save outputs
# ------------------------------------------------------------------------------

save_intermediate(list(de_table = de_table,
                       signature = signature,
                       n_pseudobulks = c(NeuralLineage = sum(pb_meta$contrast == "NeuralLineage"),
                                         Other = sum(pb_meta$contrast == "Other")),
                       presynaptic_marker_hits = present),
                  "stage2_gbm_neural_mimicry")

write_tsv(signature, here::here(PARAMS$paths$intermediate, "stage2_signature.tsv"))

snapshot_session(STAGE)
log_msg(STAGE, "Stage 2 complete.")
