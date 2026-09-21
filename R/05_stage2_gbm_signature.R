# ==============================================================================
# 05_stage2_gbm_signature.R
# Purpose: Stage 2 - characterise GBM neural-mimicry signature from Neftel et al.
# 2019 scRNA-seq data (GSE131928, Smart-seq2 platform, adult IDH-wt only).
#
# Pre-specified in Analysis Plan Section 6.2.
# Deviations documented in DEVIATIONS.md (2026-09-21 Neftel entry).
#
# Analysis design:
#   1. Load Neftel Smart-seq2 TPM matrix (5,742 adult IDH-wt cells x 20 tumors)
#   2. Reverse log-transform to linear TPM (GEO file already log-normalised)
#   3. Apply Neftel Ea > 4 filter on linear TPM (8,654 genes retained)
#   4. Compute Seurat AddModuleScore (Neftel Methods: nbin=30, ctrl=100)
#   5. Assign cellular states via Neftel Fig 3B rules:
#        - Cycling if G1S/G2M > qnorm(0.999) upper tail
#        - Non-cycling assigned to highest of 6 meta-modules
#        - Collapse MES1/MES2 -> MES, NPC1/NPC2 -> NPC
#        - Hybrids via Neftel 3-criteria (score>1, top-10%, delta>=0.3)
#   6. Pseudobulk DE: NPC-like vs (AC-like + MES-like) per tumor
#
# Contrast rationale (H2):
#   Neural-mimicry states = NPC-like + OPC-like (Neftel Fig 5A)
#   Reference: AC-like + MES-like (astrocytic/mesenchymal, non-developmental)
#   Signal captures GBM re-expression of neural developmental programs
#   that may recapitulate vulnerable RORB+ AD neuron identity loss.
#
# Outputs (results/intermediate/):
#   neftel_states.rds       Cell-level scores + state assignments
#   neftel_pseudobulk.rds   Pseudobulk counts (tumor x state)
#   neftel_deg_primary.rds  edgeR NeuralMimicry vs Other DE table
# ==============================================================================

source(here::here("R", "00_setup.R"))
STAGE <- "stage2"

# ------------------------------------------------------------------------------
# Helper 1: Load Neftel Smart-seq2 matrix + adult metadata, apply adult filter
# ------------------------------------------------------------------------------
load_neftel_adult <- function() {
  tpm_path <- here::here(PARAMS$paths$raw_data, "GSE131928",
                          "GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv.gz")
  meta_path <- here::here(PARAMS$paths$raw_data, "GSE131928",
                           "GSE131928_single_cells_tumor_name_and_adult_or_peidatric.xlsx")
  assert_that(file.exists(tpm_path), sprintf("Neftel TPM not found: %s", tpm_path))
  assert_that(file.exists(meta_path), sprintf("Neftel metadata not found: %s", meta_path))
  
  log_msg(STAGE, "Loading Neftel Smart-seq2 TPM (may take 30-60 sec)...")
  tpm_mat <- data.table::fread(tpm_path, header = TRUE, data.table = FALSE)
  log_msg(STAGE, sprintf("Raw TPM matrix: %d genes x %d cells", nrow(tpm_mat), ncol(tpm_mat) - 1))
  
  neftel_meta_full <- readxl::read_excel(meta_path, skip = 43, col_names = TRUE)
  ss2_meta <- neftel_meta_full[grepl("Smartseq2", neftel_meta_full$`processed data file`), ]
  adult_cells <- ss2_meta$`Sample name`[ss2_meta$`adult/pediatric` == "adult"]
  
  common <- intersect(colnames(tpm_mat)[-1], adult_cells)
  assert_that(length(common) > 0, "No overlap between TPM cells and adult metadata")
  log_msg(STAGE, sprintf("Adult Smart-seq2 cells retained: %d", length(common)))
  
  adult_tpm <- as.matrix(tpm_mat[, c("GENE", common)])
  rownames(adult_tpm) <- adult_tpm[, "GENE"]
  adult_tpm <- adult_tpm[, -1]
  storage.mode(adult_tpm) <- "numeric"
  
  adult_meta <- as.data.frame(ss2_meta[ss2_meta$`Sample name` %in% common, ])
  rownames(adult_meta) <- adult_meta$`Sample name`
  adult_meta <- adult_meta[colnames(adult_tpm), ]
  
  list(tpm_log = adult_tpm, meta = adult_meta)  # tpm_log = already log2(TPM/10+1)
}

# ------------------------------------------------------------------------------
# Helper 2: Reverse log transform, apply Neftel Ea filter, compute Er
# ------------------------------------------------------------------------------
preprocess_neftel <- function(tpm_log, ea_threshold = 4) {
  # GEO file is log2(TPM/10 + 1); reverse to get linear TPM
  linear_tpm <- 10 * (2^tpm_log - 1)
  Ea_linear <- log2(rowMeans(linear_tpm) + 1)
  
  genes_keep <- Ea_linear > ea_threshold
  log_msg(STAGE, sprintf("Ea > %g filter: %d/%d genes retained",
                          ea_threshold, sum(genes_keep), length(Ea_linear)))
  
  E_mat_filt <- tpm_log[genes_keep, ]
  Er_mat <- E_mat_filt - rowMeans(E_mat_filt)  # gene-centered
  list(linear_tpm = linear_tpm[genes_keep, ],
       E_filt = E_mat_filt,
       Er = Er_mat)
}

# ------------------------------------------------------------------------------
# Helper 3: AddModuleScore via Seurat (Neftel Methods replication)
# ------------------------------------------------------------------------------
score_neftel_modules <- function(linear_tpm, Er_mat, meta, modules,
                                  nbin = 30, ctrl = 100, seed = 20260101) {
  set.seed(seed)
  seu <- Seurat::CreateSeuratObject(counts = linear_tpm, meta.data = meta,
                                     min.cells = 0, min.features = 0)
  seu <- Seurat::SetAssayData(seu, assay = "RNA", layer = "data",
                               new.data = Er_mat)
  seu <- Seurat::AddModuleScore(seu, features = modules,
                                  name = paste0(names(modules), "__"),
                                  nbin = nbin, ctrl = ctrl, seed = seed)
  score_cols <- grep("__[0-9]+$", colnames(seu@meta.data), value = TRUE)
  neftel_scores <- as.matrix(seu@meta.data[, score_cols])
  colnames(neftel_scores) <- names(modules)
  neftel_scores
}

# ------------------------------------------------------------------------------
# Helper 4: Assign cellular states (Neftel Fig 3B rules)
# ------------------------------------------------------------------------------
assign_neftel_states <- function(scores, meta) {
  df <- as.data.frame(scores)
  # Rename G1/S and G2/M to R-friendly names
  colnames(df) <- gsub("/", "", colnames(df))  # G1/S -> G1S, G2/M -> G2M
  df$cell_id <- rownames(df)
  df$tumor <- meta$`tumour name`[match(df$cell_id, rownames(meta))]
  
  # Cycling: normal-fit on G1S/G2M, P<0.001 upper tail (Neftel Methods)
  fit_g1s <- MASS::fitdistr(df$G1S, "normal")$estimate
  fit_g2m <- MASS::fitdistr(df$G2M, "normal")$estimate
  thr_g1s <- fit_g1s["mean"] + qnorm(0.999) * fit_g1s["sd"]
  thr_g2m <- fit_g2m["mean"] + qnorm(0.999) * fit_g2m["sd"]
  df$cycling <- (df$G1S > thr_g1s) | (df$G2M > thr_g2m)
  
  # Assignment: highest of 6 meta-modules
  mm_cols <- c("MES2", "MES1", "AC", "OPC", "NPC1", "NPC2")
  df$top_state <- mm_cols[apply(df[, mm_cols], 1, which.max)]
  df$top_score <- apply(df[, mm_cols], 1, max)
  df$second_state <- mm_cols[apply(df[, mm_cols], 1, function(x) order(x, decreasing = TRUE)[2])]
  df$second_score <- apply(df[, mm_cols], 1, function(x) sort(x, decreasing = TRUE)[2])
  df$third_score <- apply(df[, mm_cols], 1, function(x) sort(x, decreasing = TRUE)[3])
  
  # Collapse MES1/MES2 -> MES, NPC1/NPC2 -> NPC
  collapse <- c(MES2 = "MES", MES1 = "MES", AC = "AC", OPC = "OPC",
                 NPC1 = "NPC", NPC2 = "NPC")
  df$main_state <- collapse[df$top_state]
  
  # Hybrid classification: Neftel 3 criteria
  h1 <- df$second_score > 1
  per_state_10 <- tapply(df$top_score, df$top_state, function(x) quantile(x, 0.10))
  h2 <- df$second_score > per_state_10[df$second_state]
  h3 <- (df$second_score - df$third_score) >= 0.3
  df$is_hybrid <- h1 & h2 & h3
  
  same_family <- (df$top_state %in% c("MES1","MES2") & df$second_state %in% c("MES1","MES2")) |
                  (df$top_state %in% c("NPC1","NPC2") & df$second_state %in% c("NPC1","NPC2"))
  df$is_hybrid <- df$is_hybrid & !same_family
  
  log_msg(STAGE, sprintf("Cycling: %d/%d (%.1f%%)",
                          sum(df$cycling), nrow(df), 100*mean(df$cycling)))
  log_msg(STAGE, sprintf("Hybrid:  %d/%d (%.1f%%)",
                          sum(df$is_hybrid), nrow(df), 100*mean(df$is_hybrid)))
  log_msg(STAGE, "Main state distribution:")
  print(table(df$main_state))
  
  list(df = df,
       cycling_thresholds = c(G1S_p0.001 = as.numeric(thr_g1s),
                              G2M_p0.001 = as.numeric(thr_g2m)),
       cycling_fit = list(G1S = as.list(fit_g1s), G2M = as.list(fit_g2m)))
}

# ------------------------------------------------------------------------------
# Helper 5: Pseudobulk (tumor x state) + edgeR DE for H2 contrast
# ------------------------------------------------------------------------------
compute_neftel_pseudobulk_de <- function(linear_tpm, df, min_cells = 10) {
  # Assign cells to contrast groups (H2 pre-specified)
   #   NeuralLineage = NPC + OPC (developmental mimicry)
   #   Other         = AC + MES (astrocytic/mesenchymal)
  df$contrast_group <- ifelse(df$main_state %in% c("NPC", "OPC"),
                              "NeuralLineage", "Other")
  
  # Exclude cycling and hybrid cells for cleaner signature
  keep_cells <- !df$cycling & !df$is_hybrid
  df_clean <- df[keep_cells, ]
  log_msg(STAGE, sprintf("Non-cycling non-hybrid cells for DE: %d/%d",
                          nrow(df_clean), nrow(df)))
  log_msg(STAGE, sprintf("  NeuralLineage (NPC+OPC): %d",
                          sum(df_clean$contrast_group == "NeuralLineage")))
  log_msg(STAGE, sprintf("  Other (AC+MES):          %d",
                          sum(df_clean$contrast_group == "Other")))
  
  # Pseudobulk: sum linear TPM per tumor x group
  # (Neftel dataset is TPM not counts; edgeR handles TPM-summed pseudobulk
  # via TMM normalisation - Squair et al. 2021 acceptable for smart-seq2)
  group_id <- paste(df_clean$tumor, df_clean$contrast_group, sep = "__")
  cells_in_matrix <- intersect(df_clean$cell_id, colnames(linear_tpm))
  df_clean <- df_clean[df_clean$cell_id %in% cells_in_matrix, ]
  
  # Aggregate by summing (pseudobulk)
  agg_mat <- matrix(0, nrow = nrow(linear_tpm),
                    ncol = length(unique(group_id)),
                    dimnames = list(rownames(linear_tpm),
                                    sort(unique(group_id))))
  for (g in colnames(agg_mat)) {
    cells_g <- df_clean$cell_id[group_id == g]
    if (length(cells_g) > 0) {
      agg_mat[, g] <- rowSums(linear_tpm[, cells_g, drop = FALSE])
    }
  }
  
  n_cells_per_group <- as.integer(table(group_id))
  names(n_cells_per_group) <- names(table(group_id))
  
  # Filter: keep pseudobulk samples with >= min_cells
  keep_pb <- n_cells_per_group[colnames(agg_mat)] >= min_cells
  agg_mat <- agg_mat[, keep_pb]
  
  # Sample metadata for design
  pb_meta <- data.frame(
    sample_id      = colnames(agg_mat),
    tumor          = sapply(strsplit(colnames(agg_mat), "__"), `[`, 1),
    contrast_group = factor(sapply(strsplit(colnames(agg_mat), "__"), `[`, 2),
                             levels = c("Other", "NeuralLineage")),
    n_cells        = n_cells_per_group[colnames(agg_mat)],
    stringsAsFactors = FALSE
  )
  log_msg(STAGE, sprintf("Pseudobulk samples: %d total after min_cells=%d filter",
                          ncol(agg_mat), min_cells))
  log_msg(STAGE, sprintf("  NeuralLineage: %d samples, Other: %d samples",
                          sum(pb_meta$contrast_group == "NeuralLineage"),
                          sum(pb_meta$contrast_group == "Other")))
  
  # edgeR paired design: tumor as blocking factor + contrast_group
  design <- model.matrix(~ tumor + contrast_group, data = pb_meta)
  y <- edgeR::DGEList(counts = round(agg_mat), samples = pb_meta)
  keep_g <- edgeR::filterByExpr(y, design = design)
  y <- y[keep_g, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y, method = "TMM")
  y <- edgeR::estimateDisp(y, design = design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design = design, robust = TRUE)
  
  # Contrast: NeuralLineage vs Other
  # coefficient is contrast_groupNeuralLineage (last column of design)
  qlf <- edgeR::glmQLFTest(fit, coef = "contrast_groupNeuralLineage")
  tt <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
  tt$gene <- rownames(tt)
  tt <- tt[, c("gene", "logFC", "logCPM", "F", "PValue", "FDR")]
  tt <- tt[order(tt$PValue), ]
  
  log_msg(STAGE, sprintf("Neftel DE [NeuralLineage vs Other]: %d genes, %d FDR<0.05",
                          nrow(tt), sum(tt$FDR < 0.05, na.rm = TRUE)))
  
  list(pseudobulk = agg_mat, pb_meta = pb_meta,
       de_table = tt, contrast = "NeuralLineage vs Other")
}

# ==============================================================================
# Main pipeline
# ==============================================================================

# 1. Load Neftel modules from Stage 1 preprocessing
neftel_modules <- load_intermediate("neftel_modules")
assert_that(length(neftel_modules) == 8,
             "Neftel modules should have 8 entries (MES1/2, AC, OPC, NPC1/2, G1S, G2M)")

# 2. Load + adult filter
nf <- load_neftel_adult()

# 3. Preprocess (reverse log, Ea filter, gene-center)
pp <- preprocess_neftel(nf$tpm_log)

# 4. Score with Seurat AddModuleScore
scores <- score_neftel_modules(pp$linear_tpm, pp$Er, nf$meta, neftel_modules)

# 5. Assign states
st <- assign_neftel_states(scores, nf$meta)

# 6. Save state assignments
neftel_states <- list(
  scores = scores, assignments = st$df,
  cycling_thresholds = st$cycling_thresholds,
  cycling_fit = st$cycling_fit,
  n_cells = nrow(st$df), n_tumors = length(unique(st$df$tumor)),
  filter_used = "Ea > 4 on linear TPM (Neftel default)",
  method = "Seurat AddModuleScore, nbin=30, ctrl=100, seed=20260101",
  cycling_method = "Normal-fit on G1S/G2M, P<0.001 upper tail"
)
save_intermediate(neftel_states, "neftel_states")

# 7. Pseudobulk + DE (NeuralLineage vs Other)
de_res <- compute_neftel_pseudobulk_de(pp$linear_tpm, st$df)
save_intermediate(de_res, "neftel_deg_primary")

log_msg(STAGE, "Stage 2 complete.")
