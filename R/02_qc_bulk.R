# ==============================================================================
# 02_qc_bulk.R
# Purpose: RMA normalisation, filtering, ComBat batch correction for AD bulk;
# TPM filtering for TCGA-GBM + GTEx.
# Outputs:
#   - ad_bulk_expr (batch-corrected, filtered expression matrix)
#   - ad_bulk_meta (donor metadata)
#   - gbm_bulk_expr, gbm_bulk_meta
# ==============================================================================

source(here::here("R", "00_setup.R"))
STAGE <- "qc_bulk"

# ------------------------------------------------------------------------------
# GSE48350 (HG-U133 Plus 2.0) - RMA via affy
# ------------------------------------------------------------------------------

process_gse48350 <- function() {
  log_msg(STAGE, "Processing GSE48350 (HG-U133 Plus 2.0)")
  raw_dir <- here::here(PARAMS$paths$raw_data, "GSE48350")

  # Load CEL files
  cel_files <- list.files(raw_dir, pattern = "\\.CEL(\\.gz)?$", 
                        full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  assert_that(length(cel_files) > 0, "No CEL files found for GSE48350")

  affy_batch <- affy::ReadAffy(filenames = cel_files)
  eset <- affy::rma(affy_batch)
  expr <- Biobase::exprs(eset)

  # Probe-to-gene mapping (hgu133plus2.db)
  library(hgu133plus2.db)
  probe_gene <- select(hgu133plus2.db,
                       keys = rownames(expr),
                       columns = "SYMBOL",
                       keytype = "PROBEID") %>%
    filter(!is.na(SYMBOL)) %>%
    distinct(PROBEID, .keep_all = TRUE)

  # Keep highest-mean probe per gene
  expr_df <- expr[probe_gene$PROBEID, ]
  rowmeans <- rowMeans(expr_df)
  best_probe <- probe_gene %>%
    mutate(mean_expr = rowmeans[PROBEID]) %>%
    arrange(SYMBOL, desc(mean_expr)) %>%
    distinct(SYMBOL, .keep_all = TRUE)

  gene_expr <- expr_df[best_probe$PROBEID, ]
  rownames(gene_expr) <- best_probe$SYMBOL
  
  # ---- Normalise column names to GSM IDs ----
  # RMA keeps the full CEL file name as column name (e.g. "GSM1176196_1105A-08_EC_48_Affy.CEL.gz").
  # Metadata is indexed by short GSM ID (e.g. "GSM1176196"), so we strip everything after the first "_"
  # and the .CEL(.gz) extension.
  colnames(gene_expr) <- sub("_.*", "",
                              sub("\\.CEL(\\.gz)?$", "",
                                  basename(colnames(gene_expr)),
                                  ignore.case = TRUE))

  # ---- Sample metadata (using pre-parsed GEO columns) ----
  gse <- getGEO(filename = list.files(raw_dir, pattern = "series_matrix", full.names = TRUE)[1])
  pheno <- pData(gse)

  # GSE48350 diagnosis derivation:
  # - The 80 samples with a Braak stage annotation are ALL AD patients (verified from title
  #   field where diagnosis is embedded as "_AD_"). Braak stage distinguishes disease progression.
  # - The 140 samples coded "C" in characteristics_ch1 are Caucasian controls (healthy aging).
  # - The 33 samples coded "AA" are African American controls, aged 30-48 (too young for AD
  #   analysis; excluded).
  # See DEVIATIONS.md for rationale.
  
  # Ethnicity code from characteristics_ch1 ("individual: N, X" where X in {C, AA})
  diag_codes <- sapply(strsplit(as.character(pheno$characteristics_ch1), ","), 
                       function(x) trimws(x[2]))
  
  # Derive diagnosis:
  # - Braak-annotated samples -> AD (all 80 are AD per title verification)
  # - "C"-coded samples without Braak -> Control
  # - "AA"-coded samples -> excluded (young controls, age mismatch)
  diagnosis <- rep(NA_character_, nrow(pheno))
  has_braak <- !is.na(pheno$`braak stage:ch1`)
  diagnosis[has_braak] <- "AD"
  diagnosis[!has_braak & diag_codes == "C"] <- "Control"
  # AA samples remain NA -> will be excluded
  
  meta <- data.frame(
    sample_id = rownames(pheno),
    diagnosis = diagnosis,
    age = as.numeric(pheno$`age (yrs):ch1`),
    sex = pheno$`gender:ch1`,
    brain_region = pheno$`brain region:ch1`,
    braak_stage = pheno$`braak stage:ch1`,
    apoe = pheno$`apoe genotype:ch1`,
    mmse = as.numeric(pheno$`mmse:ch1`),
    ethnicity_code = diag_codes,
    dataset = "GSE48350",
    stringsAsFactors = FALSE
  )
  
  log_msg(STAGE, sprintf("Diagnosis assignments: AD=%d, Control=%d, Excluded=%d",
                         sum(meta$diagnosis == "AD", na.rm = TRUE),
                         sum(meta$diagnosis == "Control", na.rm = TRUE),
                         sum(is.na(meta$diagnosis))))
  
  # Match sample order to expression columns
  meta <- meta[match(colnames(gene_expr), meta$sample_id), ]
  
  # ---- Filters ----
  # 1. Resolvable diagnosis
  keep_diag <- !is.na(meta$diagnosis)
  # 2. Age >= 65 (applies to both groups; AD group naturally satisfies; excludes young C samples)
  keep_age <- !is.na(meta$age) & meta$age >= 65
  # 3. Brain region: hippocampus OR entorhinal cortex
    keep_region <- grepl("hippocampus", tolower(meta$brain_region))
  # 4. For AD samples, restrict to advanced Braak stages (V, V-VI, VI)
  #    to align with the snRNA-seq analysis (Leng: Braak 0 vs Braak VI)
  advanced_braak <- meta$braak_stage %in% c("V", "V-VI", "VI")
  keep_ad_advanced <- ifelse(meta$diagnosis == "AD", advanced_braak, TRUE)
  
  keep <- keep_diag & keep_age & keep_region & keep_ad_advanced
  
  log_msg(STAGE, sprintf("Filter cascade: total=%d -> diagnosis=%d -> age65+=%d -> region=%d -> advanced Braak=%d",
                         nrow(meta), sum(keep_diag), 
                         sum(keep_diag & keep_age),
                         sum(keep_diag & keep_age & keep_region),
                         sum(keep)))
  
  meta <- meta[keep, ]
  gene_expr <- gene_expr[, keep]
  
  log_msg(STAGE, sprintf("Final GSE48350 cohort composition:"))
  log_msg(STAGE, sprintf("  Diagnosis: %s", 
                         paste(names(table(meta$diagnosis)), table(meta$diagnosis), 
                               sep = "=", collapse = ", ")))
  log_msg(STAGE, sprintf("  Region: %s",
                         paste(names(table(meta$brain_region)), table(meta$brain_region),
                               sep = "=", collapse = ", ")))
  log_msg(STAGE, sprintf("  Sex: %s",
                         paste(names(table(meta$sex)), table(meta$sex),
                               sep = "=", collapse = ", ")))

  log_msg(STAGE, sprintf("GSE48350 after QC: %d genes x %d samples", nrow(gene_expr), ncol(gene_expr)))
  list(expr = gene_expr, meta = meta)
}

# ------------------------------------------------------------------------------
# GSE36980 (HuGene 1.0 ST) - RMA via oligo
# ------------------------------------------------------------------------------

process_gse36980 <- function() {
  log_msg(STAGE, "Processing GSE36980 (HuGene 1.0 ST)")
  raw_dir <- here::here(PARAMS$paths$raw_data, "GSE36980")

  cel_files <- list.files(raw_dir, pattern = "\\.CEL(\\.gz)?$", 
                        full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  assert_that(length(cel_files) > 0, "No CEL files found for GSE36980")

  raw_data <- oligo::read.celfiles(cel_files)
  eset <- oligo::rma(raw_data)
  expr <- Biobase::exprs(eset)

  # Probe-to-gene mapping (hugene10sttranscriptcluster.db)
  library(hugene10sttranscriptcluster.db)
  probe_gene <- select(hugene10sttranscriptcluster.db,
                       keys = rownames(expr),
                       columns = "SYMBOL",
                       keytype = "PROBEID") %>%
    filter(!is.na(SYMBOL)) %>%
    distinct(PROBEID, .keep_all = TRUE)

  expr_df <- expr[probe_gene$PROBEID, ]
  rowmeans <- rowMeans(expr_df)
  best_probe <- probe_gene %>%
    mutate(mean_expr = rowmeans[PROBEID]) %>%
    arrange(SYMBOL, desc(mean_expr)) %>%
    distinct(SYMBOL, .keep_all = TRUE)

  gene_expr <- expr_df[best_probe$PROBEID, ]
  rownames(gene_expr) <- best_probe$SYMBOL
  # ---- Normalise column names to GSM IDs ----
  # oligo::rma keeps the full CEL file name as column name; metadata uses short GSM IDs.
  colnames(gene_expr) <- sub("_.*", "",
                              sub("\\.CEL(\\.gz)?$", "",
                                  basename(colnames(gene_expr)),
                                  ignore.case = TRUE))

   # ---- Sample metadata (using pre-parsed GEO columns) ----
  gse <- getGEO(filename = list.files(raw_dir, pattern = "series_matrix", full.names = TRUE)[1])
  pheno <- pData(gse)
  
  # GSE36980 diagnosis derivation:
  # Title format is "AD_XX, biological repN" or "non-AD_XX, biological repN"
  # The prefix before the first "_" indicates diagnosis.
  # Verified against pheno$title inspection during data exploration.
  title_prefix <- sub("_.*", "", pheno$title)
  diagnosis <- dplyr::case_when(
    title_prefix == "AD" ~ "AD",
    title_prefix == "non-AD" ~ "Control",
    TRUE ~ NA_character_
  )
  
  meta <- data.frame(
    sample_id = rownames(pheno),
    diagnosis = diagnosis,
    age = as.numeric(pheno$`age:ch1`),
    sex = pheno$`Sex:ch1`,
    brain_region = pheno$`tissue:ch1`,
    braak_stage = NA_character_,   # GSE36980 does not provide Braak staging
    apoe = NA_character_,
    mmse = NA_real_,
    ethnicity_code = "Japanese",   # Hokama et al. 2014 Japanese cohort
    dataset = "GSE36980",
    stringsAsFactors = FALSE
  )
  
  log_msg(STAGE, sprintf("Diagnosis assignments: AD=%d, Control=%d, Excluded=%d",
                         sum(meta$diagnosis == "AD", na.rm = TRUE),
                         sum(meta$diagnosis == "Control", na.rm = TRUE),
                         sum(is.na(meta$diagnosis))))
  
  # Match sample order to expression columns
  meta <- meta[match(colnames(gene_expr), meta$sample_id), ]
  
  # ---- Filters ----
  # 1. Resolvable diagnosis
  keep_diag <- !is.na(meta$diagnosis)
  # 2. Age >= 65 (excludes young controls; all AD samples naturally satisfy)
  keep_age <- !is.na(meta$age) & meta$age >= 65
  # 3. Brain region: hippocampus OR temporal cortex (plan Section 4.1)
  #    Frontal cortex excluded — GSE36980 has frontal cortex samples but they do not
  #    overlap with GSE48350's anatomy (hippocampus + entorhinal cortex).
  #    The only common region across cohorts is hippocampus; temporal cortex is retained
  #    as GSE36980's anatomically closest region to entorhinal cortex.
    keep_region <- grepl("hippocampus", tolower(meta$brain_region))
  
  keep <- keep_diag & keep_age & keep_region
  
  log_msg(STAGE, sprintf("Filter cascade: total=%d -> diagnosis=%d -> age65+=%d -> region=%d",
                         nrow(meta), sum(keep_diag),
                         sum(keep_diag & keep_age),
                         sum(keep)))
  
  meta <- meta[keep, ]
  gene_expr <- gene_expr[, keep]
  
  log_msg(STAGE, sprintf("Final GSE36980 cohort composition:"))
  log_msg(STAGE, sprintf("  Diagnosis: %s",
                         paste(names(table(meta$diagnosis)), table(meta$diagnosis),
                               sep = "=", collapse = ", ")))
  log_msg(STAGE, sprintf("  Region: %s",
                         paste(names(table(meta$brain_region)), table(meta$brain_region),
                               sep = "=", collapse = ", ")))
  log_msg(STAGE, sprintf("  Sex: %s",
                         paste(names(table(meta$sex)), table(meta$sex),
                               sep = "=", collapse = ", ")))

  log_msg(STAGE, sprintf("GSE36980 after QC: %d genes x %d samples", nrow(gene_expr), ncol(gene_expr)))
  list(expr = gene_expr, meta = meta)
}

# ------------------------------------------------------------------------------
# Meta-cohort assembly + ComBat batch correction
# ------------------------------------------------------------------------------

merge_ad_cohorts <- function(d1, d2) {
  log_msg(STAGE, "Merging AD cohorts at common gene symbols")
    # ---- Harmonise categorical variables across cohorts ----
  # Different GEO submissions use different conventions for the same categories.
  # Sex: GSE36980 uses "F"/"M", GSE48350 uses "female"/"male"
  # Brain region: GSE36980 uses "Hippocampus" (capital H), GSE48350 uses "hippocampus"
  harmonise_meta <- function(m) {
    m$sex <- dplyr::case_when(
      tolower(m$sex) %in% c("f", "female") ~ "female",
      tolower(m$sex) %in% c("m", "male") ~ "male",
      TRUE ~ NA_character_
    )
    m$brain_region <- tolower(m$brain_region)
    m
  }
  d1$meta <- harmonise_meta(d1$meta)
  d2$meta <- harmonise_meta(d2$meta)
  common <- intersect(rownames(d1$expr), rownames(d2$expr))
  log_msg(STAGE, sprintf("Common genes: %d", length(common)))

  expr_merged <- cbind(d1$expr[common, ], d2$expr[common, ])
  meta_merged <- rbind(d1$meta, d2$meta)
  assert_that(all(colnames(expr_merged) == meta_merged$sample_id),
              "Sample order mismatch after merge")

  # ---- Gene-level expression filter ----
  frac_expressed <- rowSums(expr_merged > PARAMS$bulk_qc$min_log2_expression) / ncol(expr_merged)
  keep <- frac_expressed >= PARAMS$bulk_qc$min_frac_samples
  expr_merged <- expr_merged[keep, ]
  log_msg(STAGE, sprintf("Genes after expression filter: %d", nrow(expr_merged)))

  # ---- Pre-ComBat PCA diagnostic ----
  pca_pre <- prcomp(t(expr_merged), scale. = TRUE)
  r_pc1_dataset_pre <- suppressWarnings(cor(pca_pre$x[, 1], as.numeric(factor(meta_merged$dataset))))
  log_msg(STAGE, sprintf("Pre-ComBat PC1 vs dataset r = %.3f", r_pc1_dataset_pre))

  # ---- ComBat: preserve diagnosis + brain region, remove dataset effect ----
    # Build model matrix dynamically: include only covariates with >= 2 levels
  # (brain_region is now uniform "hippocampus" after Yol 1 filtering to common region)
  preserve_vars <- c("diagnosis", "brain_region", "sex")
  usable_vars <- preserve_vars[sapply(preserve_vars, function(v) {
    length(unique(na.omit(meta_merged[[v]]))) >= 2
  })]
  log_msg(STAGE, sprintf("ComBat preserved covariates: %s", paste(usable_vars, collapse = " + ")))
  
  formula_str <- paste("~", paste(usable_vars, collapse = " + "))
  mod <- model.matrix(as.formula(formula_str), data = meta_merged)
  expr_combat <- sva::ComBat(dat = expr_merged,
                              batch = meta_merged$dataset,
                              mod = mod,
                              par.prior = TRUE)

  # ---- Post-ComBat PCA acceptance check ----
  pca_post <- prcomp(t(expr_combat), scale. = TRUE)
  r_pc1_dataset_post <- suppressWarnings(cor(pca_post$x[, 1], as.numeric(factor(meta_merged$dataset))))
  r_pc1_diag_post <- suppressWarnings(cor(pca_post$x[, 1], as.numeric(factor(meta_merged$diagnosis))))
  r_pc2_diag_post <- suppressWarnings(cor(pca_post$x[, 2], as.numeric(factor(meta_merged$diagnosis))))

  log_msg(STAGE, sprintf("Post-ComBat PC1 vs dataset r = %.3f", r_pc1_dataset_post))
  log_msg(STAGE, sprintf("Post-ComBat PC1 vs diagnosis r = %.3f", r_pc1_diag_post))
  log_msg(STAGE, sprintf("Post-ComBat PC2 vs diagnosis r = %.3f", r_pc2_diag_post))

  # ---- Pre-specified acceptance criterion ----
  passed <- abs(r_pc1_dataset_post) < PARAMS$combat_acceptance$max_pc1_dataset_r &&
    max(abs(r_pc1_diag_post), abs(r_pc2_diag_post)) >= PARAMS$combat_acceptance$min_pc1_or_pc2_diagnosis_r

  if (!passed) {
    log_msg(STAGE, "WARNING: ComBat did not meet pre-specified acceptance thresholds.")
    log_msg(STAGE, "This is a deviation. Record in DEVIATIONS.md before proceeding.")
  }

  list(expr = expr_combat, meta = meta_merged,
       combat_diagnostics = list(r_pc1_dataset_pre = r_pc1_dataset_pre,
                                 r_pc1_dataset_post = r_pc1_dataset_post,
                                 r_pc1_diag_post = r_pc1_diag_post,
                                 r_pc2_diag_post = r_pc2_diag_post,
                                 passed = passed))
}

# ------------------------------------------------------------------------------
# TCGA-GBM + GTEx (RNA-seq, TOIL-normalised)
# ------------------------------------------------------------------------------

process_gbm_bulk <- function() {
  log_msg(STAGE, "Processing TCGA-GBM + GTEx (Xena TOIL)")
  xena_dir <- here::here(PARAMS$paths$raw_data, "xena")

  # NOTE: file naming below is illustrative; adapt to actual Xena downloads
  gbm_file <- list.files(xena_dir, pattern = "tcga.*gbm.*tpm", ignore.case = TRUE, full.names = TRUE)[1]
  gtex_file <- list.files(xena_dir, pattern = "gtex.*tpm", ignore.case = TRUE, full.names = TRUE)[1]

  assert_that(file.exists(gbm_file), sprintf("TCGA-GBM file not found in %s", xena_dir))
  assert_that(file.exists(gtex_file), sprintf("GTEx file not found in %s", xena_dir))

  # Load (large files — consider data.table::fread)
  gbm_expr <- as.matrix(data.table::fread(gbm_file), rownames = "gene")
  gtex_expr <- as.matrix(data.table::fread(gtex_file), rownames = "gene")

  # ---- Sample metadata ----
  # TODO: separately fetch TCGA-GBM clinical (TCGAbiolinks) and GTEx phenotype
  # For now, indicative structure:
  gbm_meta <- data.frame(
    sample_id = colnames(gbm_expr),
    diagnosis = "GBM",
    brain_region = "tumor",
    sex = NA, age = NA,
    tumor_purity = NA,   # to be filled from ABSOLUTE metadata
    sample_type = "01",  # 01 = primary; will filter recurrent
    stringsAsFactors = FALSE
  )
  gbm_meta <- gbm_meta[gbm_meta$sample_type == "01", ]

  gtex_meta <- data.frame(
    sample_id = colnames(gtex_expr),
    diagnosis = "Control",
    brain_region = NA,    # TODO: hippocampus or frontal cortex BA9 from GTEx annotation
    sex = NA, age = NA,
    hardy_score = NA,
    stringsAsFactors = FALSE
  )
  # Restrict GTEx to brain hippocampus + frontal cortex BA9
  region_ok <- grepl("hippocampus|frontal_cortex_ba9|frontal cortex \\(ba9\\)",
                     tolower(gtex_meta$brain_region))
  gtex_meta <- gtex_meta[region_ok, ]
  gtex_expr <- gtex_expr[, gtex_meta$sample_id]

  # Hardy score filter
  hardy_ok <- is.na(gtex_meta$hardy_score) | gtex_meta$hardy_score <= PARAMS$bulk_qc$gtex_hardy_max
  gtex_meta <- gtex_meta[hardy_ok, ]
  gtex_expr <- gtex_expr[, gtex_meta$sample_id]

  # ---- Merge at common genes ----
  common <- intersect(rownames(gbm_expr), rownames(gtex_expr))
  expr_merged <- cbind(gbm_expr[common, gbm_meta$sample_id],
                       gtex_expr[common, gtex_meta$sample_id])
  meta_merged <- rbind(gbm_meta[, intersect(names(gbm_meta), names(gtex_meta))],
                       gtex_meta[, intersect(names(gbm_meta), names(gtex_meta))])

  # Gene filter: log2 TPM > 1 in >= 20% samples
  frac_expressed <- rowSums(expr_merged > PARAMS$bulk_qc$min_log2_tpm) / ncol(expr_merged)
  keep <- frac_expressed >= PARAMS$bulk_qc$min_frac_samples
  expr_merged <- expr_merged[keep, ]

  log_msg(STAGE, sprintf("GBM cohort after QC: %d genes x %d samples", nrow(expr_merged), ncol(expr_merged)))
  list(expr = expr_merged, meta = meta_merged)
}

# ------------------------------------------------------------------------------
# Run pipeline
# ------------------------------------------------------------------------------

d1 <- process_gse48350()
d2 <- process_gse36980()
ad_bulk <- merge_ad_cohorts(d1, d2)
save_intermediate(ad_bulk, "ad_bulk")

gbm_bulk <- process_gbm_bulk()
save_intermediate(gbm_bulk, "gbm_bulk")

snapshot_session(STAGE)
log_msg(STAGE, "Bulk QC complete.")
