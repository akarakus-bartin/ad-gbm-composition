# ==============================================================================
# 00_setup.R
# Purpose: load all packages, load parameters, set seed, define helpers.
# Sourced by every downstream script via: source("R/00_setup.R")
# ==============================================================================

suppressPackageStartupMessages({
  # Core
  library(here)
  library(yaml)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(glue)

  # I/O and data acquisition
  library(GEOquery)      # GEO downloads
  library(R.utils)       # gunzip helpers

  # Microarray
  library(affy)          # HG-U133 Plus 2.0
  library(oligo)         # HuGene 1.0 ST
  library(sva)           # ComBat

  # DE / statistics
  library(limma)
  library(edgeR)
  library(car)           # variance inflation factor

  # Single-cell
  library(Seurat)        # v5
  library(scDblFinder)   # doublet detection
  library(SingleCellExperiment)
  library(scran)
  library(scater)

  # Deconvolution
  library(BRETIGEA)
  library(MuSiC)         # reference-based deconvolution (MIND replaced due to unavailability)

  # Networks / enrichment
  library(WGCNA)
  library(fgsea)
  library(msigdbr)
  library(decoupleR)

  # Convergence tests
  library(RRHO2)

  # Plotting
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(RColorBrewer)
})

# ---- Load parameters ----
PARAMS <- yaml::read_yaml(here::here("config", "params.yml"))

# ---- Set global seed ----
set.seed(PARAMS$seed)

# ---- Ensure output directories exist ----
for (p in PARAMS$paths) {
  dir.create(here::here(p), showWarnings = FALSE, recursive = TRUE)
}

# ---- Simple structured logger ----
log_msg <- function(stage, msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] [%s] %s", ts, stage, msg)
  message(line)
  # Also append to a persistent log
  logfile <- here::here(PARAMS$paths$logs, sprintf("%s.log", stage))
  cat(line, "\n", file = logfile, append = TRUE)
}

# ---- Save/load intermediate objects with informative filenames ----
save_intermediate <- function(obj, name) {
  path <- here::here(PARAMS$paths$intermediate, paste0(name, ".rds"))
  saveRDS(obj, path)
  log_msg("io", sprintf("Saved: %s (%.1f MB)", path, file.info(path)$size / 1e6))
  invisible(path)
}

load_intermediate <- function(name) {
  path <- here::here(PARAMS$paths$intermediate, paste0(name, ".rds"))
  if (!file.exists(path)) stop(sprintf("Intermediate not found: %s", path))
  readRDS(path)
}

# ---- Assertion helper ----
assert_that <- function(condition, msg) {
  if (!isTRUE(condition)) {
    stop(sprintf("ASSERTION FAILED: %s", msg))
  }
}

# ---- Session info snapshot (called at end of each stage) ----
snapshot_session <- function(stage) {
  path <- here::here(PARAMS$paths$logs, sprintf("session_%s.txt", stage))
  sink(path); print(sessionInfo()); sink()
  log_msg(stage, sprintf("Session info written to %s", path))
}

log_msg("setup", "Environment initialised.")
log_msg("setup", sprintf("Seed: %d", PARAMS$seed))
log_msg("setup", sprintf("R version: %s", R.version.string))
