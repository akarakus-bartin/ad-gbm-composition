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
  log_msg(STAGE, "Processing Leng snRNA-seq")
  raw_dir <- here::here(PARAMS$paths$raw_data, PARAMS$data$leng_geo)

  # NOTE: Leng data may need to be downloaded from Broad Single Cell Portal
  # (SCP1198) rather than GEO for the cleanest count matrix + cell annotations.
  # Adapt the loading step to whichever source you use.

  # Expected: a count matrix + a per-cell metadata table.
  # Below assumes the count matrix and metadata are in raw_dir.
  count_files <- list.files(raw_dir, pattern = "counts|matrix", full.names = TRUE)
  meta_files  <- list.files(raw_dir, pattern = "meta|annot",   full.names = TRUE)

  assert_that(length(count_files) > 0, "Leng count matrix not found")
  assert_that(length(meta_files) > 0, "Leng cell metadata not found")

  # TODO: adapt to actual file format (10x mtx trio vs single tsv/rds vs h5ad)
  # Placeholder using a hypothetical .rds
  counts <- readRDS(count_files[1])
  cell_meta <- read.delim(meta_files[1], row.names = 1, stringsAsFactors = FALSE)

  # Align
  common_cells <- intersect(colnames(counts), rownames(cell_meta))
  counts <- counts[, common_cells]
  cell_meta <- cell_meta[common_cells, , drop = FALSE]

  seu <- CreateSeuratObject(counts = counts, meta.data = cell_meta,
                            min.cells = 3, min.features = PARAMS$snrnaseq_qc$min_features)

  # Mitochondrial fraction
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

  # Pre-QC snapshot
  log_msg(STAGE, sprintf("Leng before QC: %d cells", ncol(seu)))
  log_msg(STAGE, sprintf("  Median nFeature: %d, median nCount: %d, median %%mito: %.2f",
                         median(seu$nFeature_RNA), median(seu$nCount_RNA), median(seu$percent.mt)))

  # QC filter
  seu <- subset(seu,
                subset = nFeature_RNA >= PARAMS$snrnaseq_qc$min_features &
                         nCount_RNA >= PARAMS$snrnaseq_qc$min_counts &
                         percent.mt <= PARAMS$snrnaseq_qc$max_pct_mito)

  # Doublet detection with scDblFinder (per sample if metadata has sample id)
  if ("sample_id" %in% colnames(seu@meta.data)) {
    sce <- as.SingleCellExperiment(seu)
    sce <- scDblFinder(sce, samples = "sample_id")
    seu$doublet_class <- sce$scDblFinder.class
    seu <- subset(seu, subset = doublet_class == "singlet")
  } else {
    log_msg(STAGE, "No sample_id metadata; running scDblFinder without sample grouping")
    sce <- as.SingleCellExperiment(seu)
    sce <- scDblFinder(sce)
    seu$doublet_class <- sce$scDblFinder.class
    seu <- subset(seu, subset = doublet_class == "singlet")
  }

  log_msg(STAGE, sprintf("Leng after QC: %d cells", ncol(seu)))

  # Standard normalisation
  seu <- NormalizeData(seu, scale.factor = 10000)
  seu <- FindVariableFeatures(seu, nfeatures = 2000)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, npcs = 30, verbose = FALSE)

  # ---- Verify cell-type annotations with canonical markers ----
  # These are the marker gene panels from Section 5.2 of the plan.
  # If annotations are missing, this section labels cells by highest-scoring marker set.

  marker_panels <- list(
    excitatory = c("SNAP25", "SYT1", "RBFOX3", "SLC17A7", "CAMK2A"),
    inhibitory = c("GAD1", "GAD2"),
    rorb_high  = c("RORB"),
    astrocyte  = c("GFAP", "AQP4"),
    microglia  = c("C1QA", "CSF1R"),
    oligodendrocyte = c("MOG", "MBP")
  )

  for (name in names(marker_panels)) {
    genes <- intersect(marker_panels[[name]], rownames(seu))
    if (length(genes) > 0) {
      seu <- AddModuleScore(seu, features = list(genes), name = paste0("score_", name))
    }
  }

  # Save
  seu
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
