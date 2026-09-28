# =============================================================================
# 12_0_download_gse125583.R — rebuild the GSE125583 cohort from public data (recount3 + GEO)
#
# The cohort in results/intermediate/gse125583_bulk.rds was originally assembled interactively
# and saved by R/12_build_gse125583_bulk.R. This script rebuilds it from scratch and checks
# that it reproduces the saved object. It does NOT overwrite the locked file.
#
# Normalisation steps of the original session were not recorded, so three candidates are
# built and compared with the saved object; the one that reproduces it is reported.
#
# Output: results/intermediate/gse125583_bulk_rebuilt.rds, outputs/gse125583_rebuild_check.txt
# Needs: recount3, GEOquery, SummarizedExperiment, edgeR, limma (Bioconductor)
# =============================================================================
suppressPackageStartupMessages({
  library(recount3); library(GEOquery); library(SummarizedExperiment); library(edgeR); library(limma) })
log_lines <- character(); say <- function(...) { x <- paste0(...); message(x); log_lines <<- c(log_lines, x) }
old <- readRDS(here::here("results/intermediate/gse125583_bulk.rds"))
say("Saved cohort: ", nrow(old$expr), " genes x ", ncol(old$expr), " samples; ",
    paste(names(table(old$meta$diagnosis)), table(old$meta$diagnosis), collapse = ", "))

# ---- 1. recount3 ---------------------------------------------------------------------
proj <- available_projects(organism = "human")
p <- proj[proj$project == "SRP181886" & proj$project_type == "data_sources", ]
stopifnot(nrow(p) == 1)
rse <- create_rse(p, type = "gene")
raw <- assay(rse, "raw_counts")
gm  <- as.data.frame(rowData(rse))
cd  <- as.data.frame(colData(rse))
say("recount3 SRP181886: ", nrow(raw), " genes x ", ncol(raw), " samples (annotation ", metadata(rse)$annotation, ")")

# ---- 2. GEO phenotypes, matched by SRA experiment accession (not by position) --------
ph <- pData(getGEO("GSE125583", GSEMatrix = TRUE)[[1]])
rel <- apply(ph[, grep("^relation", names(ph)), drop = FALSE], 1, paste, collapse = " ")
ph$srx <- regmatches(rel, regexpr("SRX[0-9]+", rel))
key <- cd[["sra.experiment_acc"]]
m <- match(key, ph$srx)
say("Matched recount3 samples to GEO by SRX: ", sum(!is.na(m)), "/", length(m))
stopifnot(all(!is.na(m)))
ph <- ph[m, ]
chars <- grep(":ch1$", names(ph), value = TRUE)
for (cn in chars) { tb <- sort(table(ph[[cn]]), decreasing = TRUE)
  say("  ", cn, ": ", paste(head(names(tb), 6), head(tb, 6), sep = "=", collapse = "; ")) }

pick <- function(rx) { h <- grep(rx, chars, ignore.case = TRUE, value = TRUE); if (length(h)) h[1] else NA }
c_dx <- pick("diagnos|disease|group|status"); c_age <- pick("^age"); c_sex <- pick("sex|gender")
say("Columns used: diagnosis = ", c_dx, "; age = ", c_age, "; sex = ", c_sex)
dx  <- tolower(ph[[c_dx]])
age <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", ph[[c_age]])))
sex <- ph[[c_sex]]
c_bk <- pick("braak"); bk <- toupper(trimws(ph[[c_bk]])); is_ad <- grepl("alzheimer|^ad$", dx) & bk %in% c("V", "VI"); is_ctl <- grepl("control|normal|^ctl", dx)
sel <- (is_ad | is_ctl) & !is.na(age) & !is.na(sex) & sex != ""
say("Advanced AD = AD diagnosis and Braak V-VI. Braak in selected AD: ", paste(names(table(bk[sel & is_ad])), table(bk[sel & is_ad]), sep = "=", collapse = "; "))
say("Braak in selected controls: ", paste(names(table(bk[sel & is_ctl])), table(bk[sel & is_ctl]), sep = "=", collapse = "; "))
say("Selected: AD ", sum(sel & is_ad), ", Control ", sum(sel & is_ctl), " (saved: AD ",
    sum(old$meta$diagnosis == "AD"), ", Control ", sum(old$meta$diagnosis == "Control"), ")")
if (sum(sel & is_ad) != sum(old$meta$diagnosis == "AD") || sum(sel & is_ctl) != sum(old$meta$diagnosis == "Control"))
  say("!! Group sizes differ from the saved cohort; check the diagnosis column listing above.")

# ---- 3. candidate normalisations ------------------------------------------------------
collapse <- function(L) {
  L <- L[rowMeans(L) > 1, ]
  sym <- gm$gene_name[match(rownames(L), gm$gene_id)]; ok <- !is.na(sym); L <- L[ok, ]; sym <- sym[ok]
  o <- order(sym, -rowMeans(L)); L <- L[o, ]; sym <- sym[o]
  L <- L[!duplicated(sym), ]; rownames(L) <- sym[!duplicated(sym)]; L }
cand <- list(
  logCPM_coverage     = function() edgeR::cpm(raw, log = TRUE),
  logCPM_coverage_TMM = function() edgeR::cpm(edgeR::normLibSizes(edgeR::DGEList(raw)), log = TRUE),
  logCPM_reads        = function() edgeR::cpm(transform_counts(rse), log = TRUE))
for (pc in c(0.25, 0.5, 1, 3, 4, 5)) cand[[paste0('logCPM_cov_prior', pc)]] <- local({ q <- pc; function() edgeR::cpm(raw, log = TRUE, prior.count = q) })
cand$log2_cpm_plus1   <- function() log2(edgeR::cpm(raw) + 1)
cand$log2_cpm_plus0.5 <- function() log2(edgeR::cpm(raw) + 0.5)
compare <- function(E) {
  g <- intersect(rownames(old$expr), rownames(E))
  v <- head(g[order(-apply(old$expr[g, ], 1, var))], 2000)
  C <- cor(old$expr[v, ], E[v, ])
  best <- apply(C, 1, which.max)
  d <- max(abs(old$expr[g, ] - E[g, best]))
  list(genes_new = nrow(E), genes_common = length(g), genes_old = nrow(old$expr),
       one_to_one = !anyDuplicated(best), min_best_r = min(apply(C, 1, max)), max_abs_diff = d, best = best)
}
res <- list(); built <- list()
for (k in names(cand)) {
  E <- tryCatch(collapse(cand[[k]]()[, sel]), error = function(e) { say("  ", k, ": failed (", conditionMessage(e), ")"); NULL })
  if (is.null(E)) next
  built[[k]] <- E; res[[k]] <- compare(E)
  with(res[[k]], say(sprintf("  %-20s genes %d (common %d of %d); samples 1:1 %s; min r %.6f; max |diff| %.3g",
                             k, genes_new, genes_common, genes_old, one_to_one, min_best_r, max_abs_diff)))
}
win <- names(res)[which.min(sapply(res, `[[`, "max_abs_diff"))]
say("Closest candidate: ", win)

# ---- 4. build the rebuilt cohort object ----------------------------------------------
E <- built[[win]]
meta <- data.frame(sample_id = ph$geo_accession[sel],
                   diagnosis = factor(ifelse(is_ad[sel], "AD", "Control"), levels = c("Control", "AD")),
                   age = age[sel], sex = factor(sex[sel]), srx = ph$srx[sel], stringsAsFactors = FALSE)
colnames(E) <- meta$sample_id
rebuilt <- list(expr = E, meta = meta,
                note = paste("Rebuilt from recount3 SRP181886 + GEO GSE125583 by R/12_0; normalisation:", win))
saveRDS(rebuilt, here::here("results/intermediate/gse125583_bulk_rebuilt.rds"))

# ---- 5. functional check: repeat the A1 Model A fit on both objects ------------------
fitA <- function(obj) {
  des <- model.matrix(~ diagnosis + sex + age, data = obj$meta)
  tt <- topTable(eBayes(lmFit(obj$expr, des), trend = TRUE), coef = "diagnosisAD", number = Inf, sort.by = "none")
  tt$gene <- rownames(tt); tt }
tO <- fitA(old); tN <- fitA(rebuilt); g <- intersect(tO$gene, tN$gene)
nO <- sum(tO$adj.P.Val < 0.05); nN <- sum(tN$adj.P.Val < 0.05)
rFC <- cor(tO$logFC[match(g, tO$gene)], tN$logFC[match(g, tN$gene)])
say(sprintf("Model A (A1): DEGs saved %d, rebuilt %d; log2FC r = %.6f", nO, nN, rFC))

r <- res[[win]]
verdict <- if (r$max_abs_diff < 1e-6 && r$genes_new == r$genes_old && r$one_to_one) "IDENTICAL" else
           if (r$min_best_r > 0.999 && rFC > 0.999 && abs(nN - nO) <= 0.01 * nO) "EQUIVALENT" else "DIFFERENT"
say("VERDICT: ", verdict)
writeLines(log_lines, here::here("outputs/gse125583_rebuild_check.txt"))
message("Report: outputs/gse125583_rebuild_check.txt")
