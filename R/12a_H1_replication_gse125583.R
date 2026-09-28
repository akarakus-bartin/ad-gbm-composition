# =============================================================================
# 12a_H1_replication_gse125583.R — EXPLORATORY replication of H1 in an independent AD cohort
# PRE-DECLARED (DEVIATIONS 2026-09-28, A1) before execution:
#   AD: GSE125583 (R/12). Model A ~ diagnosis + sex + age; Model B = A + 6 BRETIGEA scores; limma-trend.
#   GBM: Model A and Model B from outputs/stage4_H1_bretigea.rds, unchanged.
#   Primary metric: shared-DEG ratio B/A (FDR < 0.05 both cohorts, same sign).
#   Labels: <= 0.20 REPLICATES; >= 0.50 DOES NOT REPLICATE; otherwise PARTIAL.
#   Secondary: AD-only DEG retention; log2FC-threshold version (|AD| > 0.5, |GBM| > 1). MuSiC not run.
# =============================================================================
suppressPackageStartupMessages({ library(limma); library(BRETIGEA) })
g  <- readRDS(here::here("results/intermediate/gse125583_bulk.rds"))
b4 <- readRDS(here::here("outputs/stage4_H1_bretigea.rds"))
bret <- as.data.frame(BRETIGEA::brainCells(inputMat = g$expr, nMarker = 50, species = "human",
                      celltypes = c("ast", "end", "mic", "neu", "oli", "opc"), scale = TRUE))
colnames(bret) <- paste0("bret_", colnames(bret))
m <- cbind(g$meta, bret[g$meta$sample_id, , drop = FALSE])
fitm <- function(f) {
  des <- model.matrix(as.formula(f), data = m); stopifnot(nrow(des) == ncol(g$expr))
  fit <- eBayes(lmFit(g$expr, des), trend = TRUE)
  tt <- topTable(fit, coef = "diagnosisAD", number = Inf, sort.by = "none"); tt$gene <- rownames(tt)
  j <- which(colnames(des) == "diagnosisAD")
  vif <- 1 / (1 - summary(lm(des[, j] ~ des[, -c(1, j)]))$r.squared)
  list(tt = tt, formula = f, vif = vif)
}
A <- fitm("~ diagnosis + sex + age")
B <- fitm(paste("~ diagnosis + sex + age +", paste(colnames(bret), collapse = " + ")))
shared <- function(a, b, lfc = FALSE) {
  gg <- intersect(a$gene, b$gene); x <- a[match(gg, a$gene), ]; y <- b[match(gg, b$gene), ]
  ok <- x$adj.P.Val < 0.05 & y$adj.P.Val < 0.05 & sign(x$logFC) == sign(y$logFC)
  if (lfc) ok <- ok & abs(x$logFC) > 0.5 & abs(y$logFC) > 1
  c(n = sum(ok), universe = length(gg))
}
lab <- function(r) if (is.na(r)) "NOT EVALUABLE" else if (r <= 0.20) "REPLICATES" else if (r >= 0.50) "DOES NOT REPLICATE" else "PARTIAL"
res <- lapply(c(fdr = FALSE, lfc = TRUE), function(l) {
  sA <- shared(A$tt, b4$models$GBM_A$tt, l); sB <- shared(B$tt, b4$models$GBM_B$tt, l)
  r <- if (sA["n"] > 0) sB["n"] / sA["n"] else NA
  list(nA = unname(sA["n"]), nB = unname(sB["n"]), universe = unname(sA["universe"]), ratio = unname(r), label = lab(r))
})
ad_ret <- sum(B$tt$adj.P.Val < 0.05) / sum(A$tt$adj.P.Val < 0.05)
A1 <- list(n = table(g$meta$diagnosis), deg_A = sum(A$tt$adj.P.Val < 0.05), deg_B = sum(B$tt$adj.P.Val < 0.05),
           ad_retention = ad_ret, vif_B = B$vif, shared = res, bretigea = bret,
           models = list(A = A, B = B), run_time = Sys.time())
saveRDS(A1, here::here("outputs/A1_H1_replication_gse125583.rds"))
message("Saved outputs/A1_H1_replication_gse125583.rds")
