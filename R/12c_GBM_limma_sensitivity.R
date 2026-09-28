# =============================================================================
# 12c_GBM_limma_sensitivity.R — GBM Stage 2 with limma-trend on log-scaled TPM, then Stage 3 (A3)
# PRE-DECLARED (DEVIATIONS 2026-09-28, A3): Neftel pseudobulk (summed TPM) scaled to per-million,
#   log2(x + 1); ~ tumour + group; eBayes(trend = TRUE, robust = TRUE); genes = those in the edgeR table.
#   Stage 3: three tests, true RRHO2 boundaries, same thresholds (DU adj p < 0.001 via x n_pixels;
#   OR >= 3 & Bonferroni x3 p < 0.001; permutation p < 0.01; seed PARAMS$seed), AD inputs plan B6 and deviation B2.
# =============================================================================
suppressPackageStartupMessages({ library(limma); library(RRHO2) })
stopifnot(exists("PARAMS"))
f <- function(p) readRDS(here::here(p))
neft <- f("results/intermediate/neftel_deg_primary.rds"); leng <- f("results/intermediate/leng_deg_primary.rds")
s1 <- f("outputs/stage1_plan_conformant.rds")
X <- neft$pseudobulk; pm <- neft$pb_meta; X <- X[, pm$sample_id]
L <- log2(sweep(X, 2, colSums(X), "/") * 1e6 + 1)
L <- L[rownames(L) %in% neft$de_table$gene, ]
des <- model.matrix(~ tumor + contrast_group, data = pm)
fit <- eBayes(lmFit(L, des), trend = TRUE, robust = TRUE)
tt <- topTable(fit, coef = "contrast_groupNeuralLineage", number = Inf, sort.by = "none")
gbm_limma <- data.frame(gene = rownames(tt), logFC = tt$logFC, PValue = tt$P.Value, FDR = tt$adj.P.Val)
univ <- intersect(leng$de_table$gene, neft$de_table$gene)
cmp <- merge(neft$de_table[, c("gene", "logFC")], gbm_limma[, c("gene", "logFC")], by = "gene")
stage3 <- function(ad, gbm) {
  P <- PARAMS$stage3; N <- P$top_n_hypergeometric
  ad <- ad[ad$gene %in% univ, ]; gbm <- gbm[gbm$gene %in% univ, ]; gg <- intersect(ad$gene, gbm$gene)
  ad <- ad[match(gg, ad$gene), ]; gbm <- gbm[match(gg, gbm$gene), ]
  rr <- RRHO2::RRHO2_initialize(data.frame(gene = gg, value = sign(ad$logFC) * -log10(ad$PValue)),
                                data.frame(gene = gg, value = sign(gbm$logFC) * -log10(gbm$PValue)),
                                stepsize = P$rrho_step_size, labels = c("AD", "GBM"), method = "hyper", log10.ind = TRUE)
  H <- rr$hypermat; br <- which(rowMeans(is.na(H)) > 0.5); bc <- which(colMeans(is.na(H)) > 0.5)
  du <- max(H[(max(br) + 1):nrow(H), 1:(min(bc) - 1)], na.rm = TRUE); du_p <- min(1, 10^(-du) * length(H))
  dd <- max(H[(max(br) + 1):nrow(H), (max(bc) + 1):ncol(H)], na.rm = TRUE)
  ad_top <- head(ad$gene[order(ad$logFC)], N); gbm_top <- head(gbm$gene[order(-gbm$logFC)], N)
  k <- length(intersect(ad_top, gbm_top)); bg <- length(univ)
  p2 <- phyper(k - 1, N, bg - N, N, lower.tail = FALSE); or <- (k * (bg - 2 * N + k)) / ((N - k)^2)
  set.seed(PARAMS$seed)
  perm <- replicate(P$permutation_n, { a <- length(intersect(sample(univ, N), sample(univ, N)))
                                        (a * (bg - 2 * N + a)) / max(1, (N - a)^2) })
  p3 <- (sum(perm >= or) + 1) / (P$permutation_n + 1)
  pass <- c(t1 = du_p < P$rrho_fdr_threshold, t2 = or >= P$hypergeometric_or_threshold && 3 * p2 < P$hypergeometric_p_threshold,
            t3 = p3 < P$permutation_p_threshold)
  nf <- sum(!pass)
  list(DU_max = du, DU_adj_p = du_p, DD_max = dd, overlap = k, OR = or, p_bonf = min(1, 3 * p2), p_perm = p3,
       pass = pass, decision = if (nf == 0) "SUPPORTED" else if (nf == 1) "INCONCLUSIVE" else "REJECTED")
}
A3 <- list(logFC_r_edgeR_vs_limma = cor(cmp$logFC.x, cmp$logFC.y), n_FDR05_limma = sum(gbm_limma$FDR < 0.05),
           n_FDR05_edgeR = sum(neft$de_table$FDR < 0.05),
           plan_B6 = stage3(s1$primary_B6$tt, gbm_limma), deviation_B2 = stage3(s1$primary_B2$tt, gbm_limma),
           gbm_limma = gbm_limma, run_time = Sys.time())
saveRDS(A3, here::here("outputs/A3_GBM_limma_sensitivity.rds"))
message("Saved outputs/A3_GBM_limma_sensitivity.rds")
