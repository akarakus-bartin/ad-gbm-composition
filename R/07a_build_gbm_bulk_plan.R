# =============================================================================
# 07a_build_gbm_bulk_plan.R — plan Table 1 / 4.1 / 5.1 GBM bulk cohort (Xena TOIL)
# Pre-declared: see DEVIATIONS 2026-09-27 Stage 4 pre-declaration (D1-D3).
# =============================================================================
suppressPackageStartupMessages(library(data.table))
stopifnot(exists("gene_meta"))
xdir     <- here::here("data/raw/xena")
ph_file  <- file.path(xdir, "TcgaTargetGTEX_phenotype.txt.gz")
tpm_file <- file.path(xdir, "TcgaTargetGtex_rsem_gene_tpm.gz")

ph <- fread(cmd = paste("gzip -dc", shQuote(ph_file)))
setnames(ph, make.names(colnames(ph)))
cat("Fenotip sütunları:", paste(colnames(ph), collapse = ", "), "\n")
study_col <- grep("study", colnames(ph), value = TRUE)[1]
type_col  <- grep("sample_type", colnames(ph), value = TRUE)[1]
dis_col   <- grep("primary.disease", colnames(ph), value = TRUE)[1]
sex_col   <- grep("gender", colnames(ph), value = TRUE)[1]

sel_tcga <- ph[[study_col]] == "TCGA" & ph[[dis_col]] == "Glioblastoma Multiforme" &
            ph[[type_col]] == "Primary Tumor"
sel_gtex <- ph[[study_col]] == "GTEX" &
            ph$detailed_category %in% c("Brain - Hippocampus", "Brain - Frontal Cortex (Ba9)")
cat("Fenotipte TCGA-GBM birincil:", sum(sel_tcga), "| GTEx hip + BA9:", sum(sel_gtex), "\n")

hdr  <- strsplit(readLines(gzcon(file(tpm_file, "rb")), n = 1), "\t")[[1]]
keep <- intersect(hdr, ph$sample[sel_tcga | sel_gtex])
cat("TPM dosyasında bulunan:", length(keep), "\n")
tpm  <- fread(cmd = paste("gzip -dc", shQuote(tpm_file)), select = c(hdr[1], keep))
mat  <- as.matrix(tpm[, -1, with = FALSE]); rownames(mat) <- sub("\\..*$", "", tpm[[1]])

# D3: Ensembl -> symbol, duplicates by highest mean
map <- setNames(gene_meta$gene_name, sub("\\..*$", "", gene_meta$gene_id))
sym <- map[rownames(mat)]
mat <- mat[!is.na(sym), ]; sym <- sym[!is.na(sym)]
o   <- order(sym, -rowMeans(mat)); mat <- mat[o, ]; sym <- sym[o]
mat <- mat[!duplicated(sym), ]; rownames(mat) <- sym[!duplicated(sym)]

# Plan 5.1 filter: exclude genes with log2 TPM >= 1 in fewer than 20% of samples
mat <- mat[rowMeans(mat >= 1) >= 0.20, ]
cat("Filtre sonrası gen:", nrow(mat), "\n")

phk  <- ph[match(colnames(mat), ph$sample)]
meta <- data.frame(
  sample_id    = colnames(mat),
  diagnosis    = factor(ifelse(phk[[study_col]] == "TCGA", "GBM", "Control"),
                        levels = c("Control", "GBM")),
  sex          = factor(phk[[sex_col]]),
  brain_region = ifelse(phk[[study_col]] == "TCGA", "tumour", phk$detailed_category),
  stringsAsFactors = FALSE
)
print(table(meta$diagnosis, meta$brain_region))
print(table(meta$diagnosis, meta$sex, useNA = "ifany"))

gbm_bulk <- list(expr = mat, meta = meta,
                 note = "Xena TOIL log2(TPM+0.001); plan cohort; built by 07a")
saveRDS(gbm_bulk, here::here("results/intermediate/gbm_bulk.rds"))
message("Saved results/intermediate/gbm_bulk.rds")

