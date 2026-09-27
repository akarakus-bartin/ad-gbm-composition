# =============================================================================
# 07b_stage4_H1_bretigea.R — H1 test, Models A and B (BRETIGEA). Plan 2.2 / 6.4.
# Pre-declared: DEVIATIONS 2026-09-27 Stage 4 D1-D12. Method 2 (Model C) NOT in this script.
# =============================================================================
suppressPackageStartupMessages({ library(limma); library(BRETIGEA) })
ad  <- readRDS(here::here("results/intermediate/ad_bulk.rds"))
gbm <- readRDS(here::here("results/intermediate/gbm_bulk.rds"))

# D11
keep_g <- !is.na(gbm$meta$sex) & gbm$meta$sex != ""
gbm$meta <- droplevels(gbm$meta[keep_g, ]); gbm$expr <- gbm$expr[, gbm$meta$sample_id]
ad$meta$diagnosis <- relevel(factor(ad$meta$diagnosis), ref = "Control")
stopifnot(identical(colnames(ad$expr), as.character(ad$meta$sample_id)) ||
          all(colnames(ad$expr) == rownames(ad$meta)))

# D9: BRETIGEA
ct <- c("ast", "end", "mic", "neu", "oli", "opc")
run_bretigea <- function(expr) {
  s <- BRETIGEA::brainCells(inputMat = expr, nMarker = 50, species = "human",
                            celltypes = ct, scale = TRUE)
  s <- as.data.frame(s); colnames(s) <- paste0("bret_", colnames(s)); s
}
ad_bret  <- run_bretigea(ad$expr)
gbm_bret <- run_bretigea(gbm$expr)

vif_diag <- function(des) {
  j <- grep("^diagnosis", colnames(des))[1]
  X <- des[, -c(1, j), drop = FALSE]
  if (ncol(X) == 0) return(1)
  1 / (1 - summary(lm(des[, j] ~ X))$r.squared)
}
fit <- function(expr, meta, formula_str) {
  des <- model.matrix(as.formula(formula_str), data = meta)
  stopifnot(nrow(des) == ncol(expr))
  f  <- eBayes(lmFit(expr, des))
  cf <- grep("^diagnosis", colnames(des), value = TRUE)[1]
  tt <- topTable(f, coef = cf, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  list(tt = tt, formula = formula_str, vif_diag = vif_diag(des))
}
ad_m  <- cbind(ad$meta,  ad_bret[colnames(ad$expr), , drop = FALSE])
gbm_m <- cbind(gbm$meta, gbm_bret[colnames(gbm$expr), , drop = FALSE])
bret_terms <- paste(colnames(ad_bret), collapse = " + ")

# D5 / D4
models <- list(
  AD_A  = fit(ad$expr,  ad_m,  "~ diagnosis + sex + age + dataset"),
  AD_B  = fit(ad$expr,  ad_m,  paste("~ diagnosis + sex + age + dataset +", bret_terms)),
  GBM_A = fit(gbm$expr, gbm_m, "~ diagnosis + sex"),
  GBM_B = fit(gbm$expr, gbm_m, paste("~ diagnosis + sex +", bret_terms))
)
# D4 sensitivity: single-region controls
for (reg in c("Brain - Frontal Cortex (Ba9)", "Brain - Hippocampus")) {
  k  <- gbm_m$diagnosis == "GBM" | gbm_m$brain_region == reg
  tag <- ifelse(grepl("Ba9", reg), "BA9", "HIP")
  models[[paste0("GBM_A_", tag)]] <- fit(gbm$expr[, k], droplevels(gbm_m[k, ]), "~ diagnosis + sex")
  models[[paste0("GBM_B_", tag)]] <- fit(gbm$expr[, k], droplevels(gbm_m[k, ]),
                                         paste("~ diagnosis + sex +", bret_terms))
}

# D6: shared DEG counts
shared <- function(a, b, lfc = FALSE) {
  g <- intersect(a$tt$gene, b$tt$gene)
  x <- a$tt[match(g, a$tt$gene), ]; y <- b$tt[match(g, b$tt$gene), ]
  ok <- x$adj.P.Val < 0.05 & y$adj.P.Val < 0.05 & sign(x$logFC) == sign(y$logFC)
  if (lfc) ok <- ok & abs(x$logFC) > 0.5 & abs(y$logFC) > 1
  list(n = sum(ok), genes = g[ok], universe = length(g))
}
decide <- function(nA, nB) {
  if (nA == 0) return(list(ratio = NA, decision = "NOT EVALUABLE (D8)"))
  r <- nB / nA
  list(ratio = r, decision = if (r >= 0.50) "REJECTED (ratio_B >= 0.50; D7)"
                             else "DEFERRED: method 2 required (D7)")
}
res <- list()
for (suf in c("", "_BA9", "_HIP")) {
  for (lfc in c(FALSE, TRUE)) {
    sA <- shared(models$AD_A, models[[paste0("GBM_A", suf)]], lfc)
    sB <- shared(models$AD_B, models[[paste0("GBM_B", suf)]], lfc)
    key <- paste0(ifelse(suf == "", "ALL", sub("_", "", suf)), ifelse(lfc, "_lfc", "_fdr"))
    res[[key]] <- c(list(nA = sA$n, nB = sB$n, universe = sA$universe,
                         genesA = sA$genes, genesB = sB$genes), decide(sA$n, sB$n))
  }
}
stage4_bret <- list(models = models, bretigea = list(ad = ad_bret, gbm = gbm_bret),
                    shared = res, primary_key = "ALL_fdr", run_time = Sys.time())
saveRDS(stage4_bret, here::here("outputs/stage4_H1_bretigea.rds"))
message("Saved outputs/stage4_H1_bretigea.rds")

