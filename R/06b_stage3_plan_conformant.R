# =============================================================================
# 06b_stage3_plan_conformant.R
# PRE-DECLARED 2026-09-27 BEFORE EXECUTION (see chat log + DEVIATIONS):
#  - Only the AD input changes. GBM input = original gbm_de. Universe = merged_de$gene.
#  - Test 1 (RRHO2): same RRHO2_initialize call as 06. Plan-conformant decision uses the max
#    -log10(P) restricted to the DU quadrant (AD bottom half x GBM top half, assuming rows = list1,
#    lists sorted descending). Adjustment x n_pixels as in 06. Global max and UD max also reported;
#    if DU and UD disagree on pass/fail, orientation must be verified before locking the decision.
#  - Test 2: identical to 06 (top-200 by logFC, Bonferroni x3, OR >= 3 & P < 0.001).
#    Sensitivity only: top-200 by signed -log10(P).
#  - Test 3: identical to 06 (1000 random pairs, seed PARAMS$seed, empirical P < 0.01).
#  - Decision: plan 2.3 (all pass = supported; >= 2 fail = rejected; 1 fail = inconclusive).
#  - Validation: first run with ORIGINAL ad_de (expected global max 6.54 at row 35, col 38).
# =============================================================================
suppressPackageStartupMessages(library(RRHO2))
stopifnot(exists("PARAMS"), exists("merged_de"), exists("ad_de"), exists("gbm_de"),
          exists("stage1_plan"))

run_stage3 <- function(ad_tab, gbm_tab, universe) {
  P <- PARAMS$stage3; N <- P$top_n_hypergeometric
  ad_tab  <- ad_tab[ad_tab$gene %in% universe, ]
  gbm_tab <- gbm_tab[gbm_tab$gene %in% universe, ]

  # ---- Test 1: RRHO2 ----
  ad_r  <- data.frame(gene = ad_tab$gene,  value = sign(ad_tab$logFC)  * -log10(ad_tab$PValue))
  gbm_r <- data.frame(gene = gbm_tab$gene, value = sign(gbm_tab$logFC) * -log10(gbm_tab$PValue))
  cg <- intersect(ad_r$gene, gbm_r$gene)
  ad_r <- ad_r[match(cg, ad_r$gene), ]; gbm_r <- gbm_r[match(cg, gbm_r$gene), ]
  rr <- RRHO2::RRHO2_initialize(list1 = ad_r, list2 = gbm_r, stepsize = P$rrho_step_size,
                                labels = c("AD_signed_logP", "GBM_signed_logP"),
                                method = "hyper", log10.ind = TRUE)
  H <- rr$hypermat; nr <- nrow(H); nc <- ncol(H); npix <- nr * nc
  qmax <- function(rows, cols) {
    M <- H[rows, cols, drop = FALSE]; v <- max(M, na.rm = TRUE)
    pos <- which(M == v, arr.ind = TRUE)[1, ] + c(min(rows), min(cols)) - 1
    list(max = v, pos = unname(pos), p_adj = 10^(-v) * npix)
  }
  top_r <- 1:floor(nr / 2); bot_r <- (floor(nr / 2) + 1):nr
  top_c <- 1:floor(nc / 2); bot_c <- (floor(nc / 2) + 1):nc
  t1 <- list(dim = c(nr, nc), rr_names = names(rr),
             global = qmax(1:nr, 1:nc), DU = qmax(bot_r, top_c), UD = qmax(top_r, bot_c),
             DD = qmax(bot_r, bot_c), UU = qmax(top_r, top_c))
  t1$passed_DU <- t1$DU$p_adj < P$rrho_fdr_threshold
  t1$passed_UD <- t1$UD$p_adj < P$rrho_fdr_threshold

  # ---- Test 2: hypergeometric (identical to 06) ----
  hyper <- function(ad_top, gbm_top) {
    ov <- intersect(ad_top, gbm_top); bg <- length(universe)
    p  <- phyper(length(ov) - 1, length(ad_top), bg - length(ad_top), length(gbm_top),
                 lower.tail = FALSE)
    a <- length(ov); b <- N - a; cc <- N - a; d <- bg - a - b - cc
    or <- (a * d) / (b * cc)
    list(overlap_n = a, overlap = ov, OR = or, p_raw = p, p_bonf = p * 3,
         passed = or >= P$hypergeometric_or_threshold && p * 3 < P$hypergeometric_p_threshold)
  }
  t2 <- hyper(head(ad_tab$gene[order(ad_tab$logFC)], N),
              head(gbm_tab$gene[order(-gbm_tab$logFC)], N))
  t2_signed <- hyper(head(ad_r$gene[order(ad_r$value)], N),
                     head(gbm_r$gene[order(-gbm_r$value)], N))

  # ---- Test 3: permutation (identical to 06) ----
  set.seed(PARAMS$seed)
  perm <- replicate(P$permutation_n, {
    s1 <- sample(universe, N); s2 <- sample(universe, N)
    a <- length(intersect(s1, s2)); b <- N - a; cc <- N - a; d <- length(universe) - a - b - cc
    (a * d) / max(1, b * cc)
  })
  p_emp <- (sum(perm >= t2$OR) + 1) / (P$permutation_n + 1)
  t3 <- list(observed_or = t2$OR, perm_median = median(perm), p_empirical = p_emp,
             passed = p_emp < P$permutation_p_threshold)

  # ---- Decision (plan 2.3) ----
  decide <- function(pass1) {
    n_fail <- sum(!c(pass1, t2$passed, t3$passed))
    if (n_fail == 0) "SUPPORTED" else if (n_fail == 1) "INCONCLUSIVE" else "REJECTED"
  }
  list(test1 = t1, test2 = t2, test2_signed = t2_signed, test3 = t3,
       decision_DU = decide(t1$passed_DU), decision_UD = decide(t1$passed_UD))
}

univ <- merged_de$gene
stage3_plan <- list(
  reproduction_original = run_stage3(ad_de, gbm_de, univ),
  plan_B6 = run_stage3(stage1_plan$primary_B6$tt, gbm_de, univ),
  deviation_B2 = run_stage3(stage1_plan$primary_B2$tt, gbm_de, univ),
  run_time = Sys.time()
)
saveRDS(stage3_plan, here::here("outputs/stage3_plan_conformant.rds"))
message("Saved outputs/stage3_plan_conformant.rds")

