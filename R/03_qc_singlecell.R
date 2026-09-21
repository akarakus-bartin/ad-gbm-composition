# ==============================================================================
# 03_qc_singlecell.R
# Purpose: snRNA-seq QC (Leng, GSE147528) and scRNA-seq QC (Neftel, GSE131928).
# Outputs:
#   - leng_seurat (QC-passed Seurat object with cell-type annotations)
#   - neftel_seurat (QC-passed Seurat object with malignant state annotations)
# ==============================================================================

source(here::here("R", "00_setup.R"))
STAGE <- "qc_singlecell"

# ------------------------------------------------------------------------------
# Leng et al. 2021 (GSE147528) - snRNA-seq
# ------------------------------------------------------------------------------

process_leng <- function() {
  log_msg(STAGE, "Processing Leng snRNA-seq (scAlign-assigned from Synapse syn21788402)")
  
  # ---- Load pre-processed EC excitatory neurons (SingleCellExperiment) ----
  # Source: Synapse syn21788402 (Leng et al. 2021 curated data).
  # This file contains fine-resolution subclustering (subclusterAssignment column)
  # with the "s" naming used in Leng et al. 2021 Fig. 2 (e.g. EC:Exc.s1, EC:Exc.s2, EC:Exc.s4).
  # We use the EC-only file because SFG serves as the anatomical control, not a discovery cohort.
  sce_path <- here::here(PARAMS$paths$raw_data, "GSE147528", "synapse_processed", 
                          "sce.EC.Exc.scAlign.rds")
  assert_that(file.exists(sce_path),
              sprintf("Leng EC excitatory subclustered SCE not found at %s. Download from Synapse syn21788402/EC_excitatoryNeurons.", sce_path))
  
  sce_ec <- readRDS(sce_path)
  log_msg(STAGE, sprintf("Loaded EC excitatory SCE: %d genes x %d cells",
                         nrow(sce_ec), ncol(sce_ec)))
  
  # ---- Validate expected metadata columns ----
  required_cols <- c("SampleID", "BraakStage", "subclusterAssignment", "clusterCellType")
  missing_cols <- setdiff(required_cols, colnames(SummarizedExperiment::colData(sce_ec)))
  assert_that(length(missing_cols) == 0,
              sprintf("Missing required metadata columns: %s", paste(missing_cols, collapse = ", ")))
  
  # ---- Report cohort composition (pre-filter) ----
  log_msg(STAGE, "Cohort composition (all EC excitatory cells):")
  log_msg(STAGE, sprintf("  Braak stages: %s",
                         paste(names(table(sce_ec$BraakStage)), 
                               table(sce_ec$BraakStage), sep = "=", collapse = ", ")))
  log_msg(STAGE, sprintf("  Subclusters:  %s",
                         paste(names(table(sce_ec$subclusterAssignment)),
                               table(sce_ec$subclusterAssignment), sep = "=", collapse = ", ")))
  
  # ---- Filter to RORB+ vulnerable subclusters (Leng Fig. 2c) ----
  # EC:Exc.s1, EC:Exc.s2, EC:Exc.s4 all express RORB and CTC-340A15.2, CTC-535M15.2.
  # These are the "vulnerable" set in Leng et al. 2021.
  # Non-vulnerable subclusters (s0, s3, s5, s6, s7, s8) are kept as internal controls
  # but placed in a separate SCE for parallel DE analysis.
  vulnerable_subclusters <- PARAMS$leng$vulnerable_subclusters
  if (is.null(vulnerable_subclusters)) {
    vulnerable_subclusters <- c("EC:Exc.s1", "EC:Exc.s2", "EC:Exc.s4")
  }
  
  vuln_cells    <- sce_ec$subclusterAssignment %in% vulnerable_subclusters
  nonvuln_cells <- !vuln_cells
  
  sce_vuln    <- sce_ec[, vuln_cells]
  sce_nonvuln <- sce_ec[, nonvuln_cells]
  
  log_msg(STAGE, sprintf("Vulnerable RORB+ subclusters (s1/s2/s4): %d cells", ncol(sce_vuln)))
  log_msg(STAGE, sprintf("Non-vulnerable control subclusters:      %d cells", ncol(sce_nonvuln)))
  
  # ---- Report Braak x Sample x Subcluster in vulnerable subset ----
  log_msg(STAGE, "Vulnerable subset — cells per donor x subcluster:")
  cross <- as.data.frame.matrix(table(sce_vuln$SampleID, sce_vuln$subclusterAssignment))
  cross <- cross[rowSums(cross) > 0, , drop = FALSE]
  # Print row by row so it fits in log
  for (donor in rownames(cross)) {
    braak <- unique(as.character(sce_vuln$BraakStage[sce_vuln$SampleID == donor]))[1]
    log_msg(STAGE, sprintf("  %s (Braak %s): %s",
                           donor, braak,
                           paste(colnames(cross), cross[donor, ], sep = "=", collapse = ", ")))
  }
  
  list(sce_vulnerable = sce_vuln,
       sce_nonvulnerable = sce_nonvuln,
       vulnerable_subclusters = vulnerable_subclusters)
}

# ------------------------------------------------------------------------------
# Neftel et al. 2019 (GSE131928) - scRNA-seq
# ------------------------------------------------------------------------------

process_neftel <- function() {
  log_msg(STAGE, "Processing Neftel scRNA-seq")
  raw_dir <- here::here(PARAMS$paths$raw_data, PARAMS$data$neftel_geo)

  # NOTE: Neftel data has two batches (10x + Smart-seq2). Suggest starting with 10x
  # subset for consistency. Metadata is available at the Broad SCP (SCP503) with
  # per-cell malignant state assignments and meta-module scores.

  count_files <- list.files(raw_dir, pattern = "counts|matrix|expression",
                            full.names = TRUE, recursive = TRUE)
  meta_files  <- list.files(raw_dir, pattern = "meta|annot|module|classification",
                            full.names = TRUE, recursive = TRUE)

  assert_that(length(count_files) > 0, "Neftel count matrix not found")
  assert_that(length(meta_files) > 0, "Neftel cell metadata not found")

  # TODO: adapt to actual file format
  counts <- readRDS(count_files[1])           # or read10X, read.table depending on format
  cell_meta <- read.delim(meta_files[1], row.names = 1, stringsAsFactors = FALSE)

  common_cells <- intersect(colnames(counts), rownames(cell_meta))
  counts <- counts[, common_cells]
  cell_meta <- cell_meta[common_cells, , drop = FALSE]

  seu <- CreateSeuratObject(counts = counts, meta.data = cell_meta,
                            min.cells = 3, min.features = PARAMS$scrnaseq_qc$min_features)

  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

  log_msg(STAGE, sprintf("Neftel before QC: %d cells", ncol(seu)))

  seu <- subset(seu,
                subset = nFeature_RNA >= PARAMS$scrnaseq_qc$min_features &
                         nCount_RNA >= PARAMS$scrnaseq_qc$min_counts &
                         percent.mt <= PARAMS$scrnaseq_qc$max_pct_mito)

  log_msg(STAGE, sprintf("Neftel after QC: %d cells", ncol(seu)))

  # Restrict to malignant cells (per original annotation)
  if ("malignant_status" %in% colnames(seu@meta.data)) {
    seu <- subset(seu, subset = malignant_status == "malignant")
    log_msg(STAGE, sprintf("Malignant cells only: %d", ncol(seu)))
  } else {
    log_msg(STAGE, "WARNING: no malignant_status column; verify annotation before proceeding")
  }

  # Verify state assignments if provided
  if ("cellular_tumor_state" %in% colnames(seu@meta.data)) {
    log_msg(STAGE, "Neftel cellular states present:")
    print(table(seu$cellular_tumor_state))
  } else {
    log_msg(STAGE, "WARNING: no cellular_tumor_state column")
    log_msg(STAGE, "  You will need to re-derive Neftel meta-module scores; see Section 5.2 of plan.")
  }

  # Normalisation (standard, though Neftel is Smart-seq2 for a subset — TPM would be more appropriate)
  seu <- NormalizeData(seu, scale.factor = 10000)
  seu <- FindVariableFeatures(seu, nfeatures = 2000)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, npcs = 30, verbose = FALSE)

  seu
}

# ------------------------------------------------------------------------------
# Run pipeline
# ------------------------------------------------------------------------------

leng_seurat <- process_leng()
save_intermediate(leng_seurat, "leng_seurat")

neftel_seurat <- process_neftel()
save_intermediate(neftel_seurat, "neftel_seurat")

snapshot_session(STAGE)
log_msg(STAGE, "Single-cell QC complete.")
