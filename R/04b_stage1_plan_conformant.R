# =============================================================================
# 04b_stage1_plan_conformant.R
# Plan-conformant Stage 1 (Analysis_Plan_v1 §6.1 + DEVIATIONS 2026-09-21 + audit 2026-09-27)
#
# PRE-DECLARED BEFORE EXECUTION:
#   Primary contrast : Braak VI vs 0 (plan 6.1). Secondary: Braak II vs 0 (DEVIATIONS 09-21).
#   Primary method   : edgeR QLF on donor-level pseudobulk (s1+s2+s4 summed), ~ braak + age
#   Sensitivity      : voom + duplicateCorrelation, ~ braak + subcluster + age, block = donor
#   Signature        : FDR < 0.05 AND log2FC < -0.5 (plan 6.1)
#   Stage 3 ranking  : sign(logFC) * -log10(PValue) from PRIMARY method
#   Requires in memory: pb (counts + meta incl. age), from build_leng_pseudobulk()
#
# NOTE: first executed 2026-09-27 18:46 by pasting into the console; this file was
#       written afterwards with identical code. See DEVIATIONS.md timing note.
# =============================================================================
suppressPackageStartupMessages({ library(edgeR); library(limma) })
stopifnot(exists("pb"),
          all(c("sample_id", "donor", "braak_stage", "subcluster", "age") %in% colnames(pb$meta)))

run_primary <- function(ref, test) {
  m   <- pb$meta[pb$meta$braak_stage %in% c(ref, test), ]
  cnt <- pb$counts[, m$sample_id, drop = FALSE]
  don <- unique(m[, c("donor", "braak_stage", "age")])
  stopifnot(!any(duplicated(don$donor)))
  cd  <- sapply(don$donor, function(d) rowSums(cnt[, m$donor == d, drop = FALSE]))
  don$braak <- factor(don$braak_stage, levels = c(ref, test))
  des <- model.matrix(~ braak + age, data = don)
  y   <- DGEList(cd)
  y   <- y[filterByExpr(y, design = des), , keep.lib.sizes = FALSE]
  y   <- normLibSizes(y)
  y   <- estimateDisp(y, des, robust = TRUE)
  fit <- glmQLFit(y, des, robust = TRUE)
  tt  <- topTags(glmQLFTest(fit, coef = paste0("braak", test)), n = Inf, sort.by = "none")$table
  tt$gene   <- rownames(tt)
  tt$signed <- sign(tt$logFC) * -log10(tt$PValue)
  list(tt = tt, donors = table(don$braak), resid_df = nrow(des) - ncol(des),
       coefs = colnames(des), r_braak_age = cor(as.numeric(don$braak), don$age))
}

run_sensitivity <- function(ref, test) {
  m   <- pb$meta[pb$meta$braak_stage %in% c(ref, test), ]
  cnt <- pb$counts[, m$sample_id, drop = FALSE]
  m$braak      <- factor(m$braak_stage, levels = c(ref, test))
  m$subcluster <- factor(make.names(m$subcluster))
  des <- model.matrix(~ braak + subcluster + age, data = m)
  y   <- DGEList(cnt)
  y   <- y[filterByExpr(y, design = des), , keep.lib.sizes = FALSE]
  y   <- normLibSizes(y)
  v   <- voom(y, des)
  cf  <- duplicateCorrelation(v, des, block = m$donor)
  v   <- voom(y, des, block = m$donor, correlation = cf$consensus)
  cf  <- duplicateCorrelation(v, des, block = m$donor)
  fit <- eBayes(lmFit(v, des, block = m$donor, correlation = cf$consensus), robust = TRUE)
  tt  <- topTable(fit, coef = paste0("braak", test), number = Inf, sort.by = "none")
  tt$gene   <- rownames(tt)
  tt$signed <- sign(tt$logFC) * -log10(tt$P.Value)
  list(tt = tt, icc = cf$consensus)
}

signature <- function(tt, fdr_col) tt$gene[tt[[fdr_col]] < 0.05 & tt$logFC < -0.5]

stage1_plan <- list(
  primary_B6   = run_primary("0", "6"),
  primary_B2   = run_primary("0", "2"),
  sens_B6      = run_sensitivity("0", "6"),
  sens_B2      = run_sensitivity("0", "2"),
  declared     = "See header; declared before execution.",
  run_time     = Sys.time()
)
stage1_plan$signature_B6 <- signature(stage1_plan$primary_B6$tt, "FDR")
stage1_plan$signature_B2 <- signature(stage1_plan$primary_B2$tt, "FDR")

saveRDS(stage1_plan, here::here("outputs/stage1_plan_conformant.rds"))
message("Saved outputs/stage1_plan_conformant.rds")

