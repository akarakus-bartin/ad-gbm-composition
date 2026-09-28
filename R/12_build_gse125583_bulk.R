# =============================================================================
# 12_build_gse125583_bulk.R — save the GSE125583 bulk cohort used in the exploratory
# analyses to disk (it previously existed only in the R session; see DEVIATIONS).
# Sample selection identical to the 2026-09-27 exploratory analysis.
# =============================================================================
need <- c("log_cpm", "geo_clinical", "advanced_ad", "ctrl_samples", "gene_meta")
miss <- need[!vapply(need, exists, logical(1))]
if (length(miss)) stop("Not in memory: ", paste(miss, collapse = ", "),
                       ". The loading code for GSE125583 must be reconstructed first.")
age_col <- grep("age", colnames(geo_clinical), ignore.case = TRUE, value = TRUE)[1]
message("Age column used: ", age_col)
age <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", geo_clinical[[age_col]])))
sex <- geo_clinical$geo_sex
sel <- (advanced_ad | ctrl_samples) & !is.na(age) & !is.na(sex)
ids <- colnames(log_cpm)
if (is.null(ids)) {
  rn <- rownames(geo_clinical)
  ids <- if (!is.null(rn) && !identical(rn, as.character(seq_len(nrow(geo_clinical))))) rn else sprintf("GSE125583_%03d", seq_len(ncol(log_cpm)))
  message("log_cpm has no column names; using ", if (identical(ids, rn)) "rownames(geo_clinical)" else "positional IDs")
}
stopifnot(length(ids) == ncol(log_cpm), nrow(geo_clinical) == ncol(log_cpm))
E <- log_cpm[, sel]; colnames(E) <- ids[sel]
E <- E[rowMeans(E) > 1, ]
sym <- gene_meta$gene_name[match(rownames(E), gene_meta$gene_id)]
keep <- !is.na(sym); E <- E[keep, ]; sym <- sym[keep]
o <- order(sym, -rowMeans(E)); E <- E[o, ]; sym <- sym[o]
E <- E[!duplicated(sym), ]; rownames(E) <- sym[!duplicated(sym)]
meta <- data.frame(sample_id = colnames(E),
                   diagnosis = factor(ifelse(advanced_ad[sel], "AD", "Control"), levels = c("Control", "AD")),
                   age = age[sel], sex = factor(sex[sel]), stringsAsFactors = FALSE)
gse125583_bulk <- list(expr = E, meta = meta,
  note = "recount3 coverage -> log-CPM (in-session object log_cpm); Advanced AD vs Control; rowMeans > 1; built by R/12")
saveRDS(gse125583_bulk, here::here("results/intermediate/gse125583_bulk.rds"))
print(table(meta$diagnosis, meta$sex)); message("Genes: ", nrow(E), " | samples: ", ncol(E))
