# =============================================================================
# 12b_DD_decomposition.R — DESCRIPTIVE 2x2 decomposition of the original DD signal (A2)
# PRE-DECLARED (DEVIATIONS 2026-09-28, A2): Braak II vs 0 fixed. Cells:
#   (no donor, no age)  = original (leng_deg_primary)
#   (no donor, age)     = outputs/age_sensitivity_results.rds $leng_de_agecov
#   (donor, no age)     = donor-summed edgeR ~ braak (sensitivity_analysis_full_state $pseudoreplication$ttB)
#   (donor, age)        = stage1_plan_conformant $primary_B2
# Per cell: RRHO2 DD and DU maxima with true (sign-change) boundaries; top-200 DD overlap and BH p.
# No decision; all four cells reported.
# =============================================================================
suppressPackageStartupMessages(library(RRHO2))
f <- function(p) readRDS(here::here(p))
leng <- f("results/intermediate/leng_deg_primary.rds"); neft <- f("results/intermediate/neftel_deg_primary.rds")
age  <- f("outputs/age_sensitivity_results.rds");      st   <- f("outputs/sensitivity_analysis_full_state.rds")
s1   <- f("outputs/stage1_plan_conformant.rds")
cells <- list(`no donor, no age` = leng$de_table, `no donor, age` = age$leng_de_agecov$de_table,
              `donor, no age` = st$pseudoreplication$ttB, `donor, age` = s1$primary_B2$tt)
univ <- intersect(leng$de_table$gene, neft$de_table$gene); message("Universe: ", length(univ))
sc <- function(t) setNames(sign(t$logFC) * -log10(t$PValue), t$gene)
gb <- sc(neft$de_table)
one <- function(t) {
  a <- sc(t); gg <- intersect(intersect(names(a), names(gb)), univ)
  rr <- RRHO2::RRHO2_initialize(data.frame(gene = gg, value = a[gg]), data.frame(gene = gg, value = gb[gg]),
                                stepsize = 100, labels = c("AD", "GBM"), method = "hyper", log10.ind = TRUE)
  H <- rr$hypermat; br <- which(rowMeans(is.na(H)) > 0.5); bc <- which(colMeans(is.na(H)) > 0.5)
  dn_r <- (max(br) + 1):nrow(H); up_c <- 1:(min(bc) - 1); dn_c <- (max(bc) + 1):ncol(H)
  N <- length(gg); ad_dn <- head(gg[order(a[gg])], 200); gb_dn <- head(gg[order(gb[gg])], 200)
  k <- length(intersect(ad_dn, gb_dn))
  c(DD_max = max(H[dn_r, dn_c], na.rm = TRUE), DU_max = max(H[dn_r, up_c], na.rm = TRUE),
    DD_top200_overlap = k, DD_top200_expected = 200 * 200 / N,
    DD_top200_p = phyper(k - 1, 200, N - 200, 200, lower.tail = FALSE), genes = N)
}
A2 <- do.call(rbind, lapply(cells, one))
saveRDS(list(table = A2, run_time = Sys.time()), here::here("outputs/A2_DD_decomposition.rds"))
message("Saved outputs/A2_DD_decomposition.rds")
