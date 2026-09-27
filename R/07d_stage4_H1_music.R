# =============================================================================
# 07d_stage4_H1_music.R — H1 method 2 (MuSiC), Model C (+ exploratory D) and final H1 decision
# Pre-declared: DEVIATIONS 2026-09-27 D13-D24. Written and committed before matrix_selected.csv
# was inspected. Model A and BRETIGEA taken from outputs/stage4_H1_bretigea.rds (07b).
# =============================================================================
suppressPackageStartupMessages({ library(limma); library(MuSiC); library(SingleCellExperiment); library(Matrix) })

# ---- Reference (D13-D15, D20-D21) ----
sel <- data.table::fread(here::here("outputs/allen_m1_reference_cells_selected.csv"))
dt  <- data.table::fread(here::here("data/raw/allen_m1/matrix_selected.csv"))
ids <- dt[[1]]; stopifnot(setequal(ids, sel$sample_name))
cnt <- t(as.matrix(dt[, -1, with = FALSE])); colnames(cnt) <- ids; rm(dt); gc()
cnt <- Matrix(cnt, sparse = TRUE)
sel <- sel[match(colnames(cnt), sel$sample_name), ]
sce <- SingleCellExperiment(assays = list(counts = cnt),
                            colData = data.frame(ref_class = sel$ref_class,
                                                 donor = sel$external_donor_name_label))
cat("Referans:", nrow(sce), "gen x", ncol(sce), "hücre\n"); print(table(sce$ref_class, sce$donor))

# ---- Bulk cohorts: identical to 07b (D11) ----
ad  <- readRDS(here::here("results/intermediate/ad_bulk.rds"))
gbm <- readRDS(here::here("results/intermediate/gbm_bulk.rds"))
keep_g <- !is.na(gbm$meta$sex) & gbm$meta$sex != ""
gbm$meta <- droplevels(gbm$meta[keep_g, ]); gbm$expr <- gbm$expr[, gbm$meta$sample_id]
ad$meta$diagnosis <- relevel(factor(ad$meta$diagnosis), ref = "Control")
b4 <- readRDS(here::here("outputs/stage4_H1_bretigea.rds"))

# ---- MuSiC (D16, D23) ----
run_music <- function(expr_log, offset) {
  lin <- pmax(2^expr_log - offset, 0)
  g   <- intersect(rownames(lin), rownames(sce))
  cat("  ortak gen:", length(g), "\n")
  r <- music_prop(bulk.mtx = lin[g, ], sc.sce = sce[g, ], clusters = "ref_class",
                  samples = "donor", select.ct = c("neu", "ast", "mic", "oli", "opc", "end"),
                  verbose = FALSE)
  p <- as.data.frame(r$Est.prop.weighted); colnames(p) <- paste0("mus_", colnames(p)); p
}
cat("MuSiC AD:\n");  ad_mus  <- run_music(ad$expr, 0)
cat("MuSiC GBM:\n"); gbm_mus <- run_music(gbm$expr, 0.001)

# ---- Model C (D17, D22) and exploratory D (D24) ----
fit <- function(expr, meta, formula_str) {
  des <- model.matrix(as.formula(formula_str), data = meta); stopifnot(nrow(des) == ncol(expr))
  f  <- eBayes(lmFit(expr, des)); cf <- grep("^diagnosis", colnames(des), value = TRUE)[1]
  tt <- topTable(f, coef = cf, number = Inf, sort.by = "none"); tt$gene <- rownames(tt)
  list(tt = tt, formula = formula_str)
}
mus_terms <- function(p) {
  v <- setdiff(colnames(p), "mus_end")                       # D17
  v[sapply(p[, v, drop = FALSE], function(x) var(x) > 0)]    # D22
}
ad_m  <- cbind(ad$meta,  ad_mus[colnames(ad$expr), , drop = FALSE],  b4$bretigea$ad[colnames(ad$expr), , drop = FALSE])
gbm_m <- cbind(gbm$meta, gbm_mus[colnames(gbm$expr), , drop = FALSE], b4$bretigea$gbm[colnames(gbm$expr), , drop = FALSE])
tA <- paste(mus_terms(ad_mus), collapse = " + "); tG <- paste(mus_terms(gbm_mus), collapse = " + ")
bret <- paste(colnames(b4$bretigea$ad), collapse = " + ")

m <- b4$models
m$AD_C  <- fit(ad$expr,  ad_m,  paste("~ diagnosis + sex + age + dataset +", tA))
m$GBM_C <- fit(gbm$expr, gbm_m, paste("~ diagnosis + sex +", tG))
m$AD_D  <- fit(ad$expr,  ad_m,  paste("~ diagnosis + sex + age + dataset +", bret, "+", tA))
m$GBM_D <- fit(gbm$expr, gbm_m, paste("~ diagnosis + sex +", bret, "+", tG))
for (reg in c("Brain - Frontal Cortex (Ba9)", "Brain - Hippocampus")) {
  k <- gbm_m$diagnosis == "GBM" | gbm_m$brain_region == reg; tag <- ifelse(grepl("Ba9", reg), "BA9", "HIP")
  m[[paste0("GBM_C_", tag)]] <- fit(gbm$expr[, k], droplevels(gbm_m[k, ]), paste("~ diagnosis + sex +", tG))
}

# ---- Shared DEGs (D6) and final H1 decision (plan 2.2) ----
shared <- function(a, b, lfc = FALSE) {
  g <- intersect(a$tt$gene, b$tt$gene); x <- a$tt[match(g, a$tt$gene), ]; y <- b$tt[match(g, b$tt$gene), ]
  ok <- x$adj.P.Val < 0.05 & y$adj.P.Val < 0.05 & sign(x$logFC) == sign(y$logFC)
  if (lfc) ok <- ok & abs(x$logFC) > 0.5 & abs(y$logFC) > 1
  sum(ok)
}
h1 <- function(rB, rC) {
  if (is.na(rB) || is.na(rC)) return("NOT EVALUABLE")
  if (rB <= 0.20 && rC <= 0.20) "SUPPORTED" else if (rB >= 0.50 || rC >= 0.50) "REJECTED" else "INCONCLUSIVE"
}
res <- list()
for (suf in c("", "_BA9", "_HIP")) for (lfc in c(FALSE, TRUE)) {
  nA <- shared(m$AD_A, m[[paste0("GBM_A", suf)]], lfc)
  nB <- shared(m$AD_B, m[[paste0("GBM_B", suf)]], lfc)
  nC <- shared(m$AD_C, m[[paste0("GBM_C", suf)]], lfc)
  rB <- if (nA > 0) nB / nA else NA; rC <- if (nA > 0) nC / nA else NA
  key <- paste0(ifelse(suf == "", "ALL", sub("_", "", suf)), ifelse(lfc, "_lfc", "_fdr"))
  res[[key]] <- list(nA = nA, nB = nB, nC = nC, ratio_B = rB, ratio_C = rC, decision = h1(rB, rC))
}
res$exploratory_D <- list(nD = shared(m$AD_D, m$GBM_D))

stage4_music <- list(proportions = list(ad = ad_mus, gbm = gbm_mus), models = m[c("AD_C", "GBM_C", "AD_D", "GBM_D",
                     "GBM_C_BA9", "GBM_C_HIP")], shared = res, primary_key = "ALL_fdr", run_time = Sys.time())
saveRDS(stage4_music, here::here("outputs/stage4_H1_music.rds"))
message("Saved outputs/stage4_H1_music.rds")

