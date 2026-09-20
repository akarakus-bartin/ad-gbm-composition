# ==============================================================================
# 01_download_data.R
# Purpose: acquire all raw datasets (GEO + UCSC Xena TOIL).
# Run once, then downstream scripts read from data/raw/.
# ==============================================================================

source(here::here("R", "00_setup.R"))
STAGE <- "download"

raw_dir <- here::here(PARAMS$paths$raw_data)

# ------------------------------------------------------------------------------
# GEO downloads
# ------------------------------------------------------------------------------

download_geo <- function(gse_id) {
  target <- file.path(raw_dir, gse_id)
  if (dir.exists(target) && length(list.files(target)) > 0) {
    log_msg(STAGE, sprintf("Skipping %s (already present).", gse_id))
    return(target)
  }
  dir.create(target, showWarnings = FALSE, recursive = TRUE)
  log_msg(STAGE, sprintf("Downloading %s ...", gse_id))

  # Series matrix (metadata) + raw supplementary files
  getGEO(gse_id, destdir = target, GSEMatrix = TRUE, getGPL = FALSE)
  getGEOSuppFiles(gse_id, baseDir = raw_dir, makeDirectory = FALSE, fetch_files = TRUE)

  log_msg(STAGE, sprintf("Done: %s", gse_id))
  target
}

for (gse in c(PARAMS$data$leng_geo,
              PARAMS$data$neftel_geo,
              PARAMS$data$ad_bulk_1_geo,
              PARAMS$data$ad_bulk_2_geo)) {
  tryCatch(
    download_geo(gse),
    error = function(e) log_msg(STAGE, sprintf("FAILED %s: %s", gse, conditionMessage(e)))
  )
}

# ------------------------------------------------------------------------------
# UCSC Xena TOIL (TCGA-GBM + GTEx harmonised)
# ------------------------------------------------------------------------------
# NOTE: manual download strongly recommended for large files.
# Below is a wrapper; if it fails, download the URLs directly with wget/curl
# and place under data/raw/xena/.

xena_dir <- file.path(raw_dir, "xena")
dir.create(xena_dir, showWarnings = FALSE, recursive = TRUE)

download_xena <- function(url, out_name) {
  out_gz <- file.path(xena_dir, paste0(out_name, ".gz"))
  out_txt <- file.path(xena_dir, out_name)
  if (file.exists(out_txt)) {
    log_msg(STAGE, sprintf("Skipping Xena %s (already present).", out_name))
    return(out_txt)
  }
  log_msg(STAGE, sprintf("Downloading Xena %s ...", url))
  download.file(url, destfile = out_gz, method = "libcurl", mode = "wb")
  R.utils::gunzip(out_gz, destname = out_txt, remove = TRUE)
  log_msg(STAGE, sprintf("Extracted to %s", out_txt))
  out_txt
}

# TCGA-GBM + GTEx unified expression matrix (large file, several GB)
# NOTE: Consult the current UCSC Xena URLs before running; the ones in params.yml
# may need updating. Alternative: use the smaller cohort-specific downloads.
tryCatch({
  download_xena(PARAMS$data$xena_gbm_url, "tcga_gbm_tpm.tsv")
}, error = function(e) log_msg(STAGE, sprintf("Xena TCGA-GBM FAILED: %s (download manually)", conditionMessage(e))))

tryCatch({
  download_xena(PARAMS$data$xena_gtex_url, "gtex_tpm.tsv")
}, error = function(e) log_msg(STAGE, sprintf("Xena GTEx FAILED: %s (download manually)", conditionMessage(e))))

# ------------------------------------------------------------------------------
# TCGA-GBM clinical metadata (for tumor purity filtering in Stage 4)
# ------------------------------------------------------------------------------
# Alternative: use TCGAbiolinks::GDCquery_clinic

log_msg(STAGE, "Downloads complete. Verify file listing below:")
for (d in list.dirs(raw_dir)) {
  files <- list.files(d, full.names = TRUE)
  if (length(files) > 0) {
    log_msg(STAGE, sprintf("  %s: %d files, %.1f MB total",
                           d, length(files), sum(file.info(files)$size, na.rm = TRUE) / 1e6))
  }
}

snapshot_session(STAGE)
