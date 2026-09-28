# 10_facts_scenarioB.R — single source of numbers for the Scenario B manuscript
f <- function(p) readRDS(here::here(p))
s1 <- f("outputs/stage1_plan_conformant.rds"); s3 <- f("outputs/stage3_plan_conformant.rds")
b4 <- f("outputs/stage4_H1_bretigea.rds");     m4 <- f("outputs/stage4_H1_music.rds")
nc <- f("outputs/stage4_negative_control.rds"); st <- f("outputs/sensitivity_analysis_full_state.rds")
ad <- f("results/intermediate/ad_bulk.rds");   gb <- f("results/intermediate/gbm_bulk.rds")
gm <- gb$meta[!is.na(gb$meta$sex) & gb$meta$sex != "", ]
cap  <- function(p) signif(min(1, p), 3)
ndeg <- function(mm, k) sum(mm[[k]]$tt$adj.P.Val < 0.05)
t3 <- function(x, lab) sprintf(
  "| %s | DU %.2f (adj p %s) | DD %.2f (adj p %s) | OR %.2f, Bonf p %s | emp p %.3f | **%s** |", lab,
  x$test1$DU$max, cap(x$test1$DU$p_adj), x$test1$DD$max, cap(x$test1$DD$p_adj),
  x$test2$OR, cap(x$test2$p_bonf), x$test3$p_empirical, x$decision_DU)
tab <- function(v) paste(names(table(v)), table(v), sep = " ", collapse = ", ")
sh <- m4$shared$ALL_fdr; dd <- st$validation_composition$dd_table; S <- dd$Category == "Stress"
facts <- c(
"# Olgu sayfası — Senaryo B (R/10_facts_scenarioB.R ile üretildi; elle düzenlemeyin)",
sprintf("Üretim: %s | git HEAD: %s", format(Sys.time(), "%Y-%m-%d %H:%M"), system("git rev-parse --short HEAD", intern = TRUE)),
"", "## H1 (birincil; Stage 4)",
sprintf("- AD kohortu: n = %d (%s); bölge: %s; veri seti: %s", nrow(ad$meta), tab(ad$meta$diagnosis),
        tab(ad$meta$brain_region), tab(ad$meta$dataset)),
sprintf("- GBM kohortu: n = %d (%s); kontrol bölgesi: %s", nrow(gm), tab(gm$diagnosis),
        tab(gm$brain_region[gm$diagnosis == "Control"])),
sprintf("- DEG FDR<0,05 — AD: A %d, B %d, C %d | GBM: A %d, B %d, C %d",
        ndeg(b4$models, "AD_A"), ndeg(b4$models, "AD_B"), ndeg(m4$models, "AD_C"),
        ndeg(b4$models, "GBM_A"), ndeg(b4$models, "GBM_B"), ndeg(m4$models, "GBM_C")),
sprintf("- Paylaşılan DEG: A %d, B %d, C %d; oran_B %.3f, oran_C %.3f -> **%s**",
        sh$nA, sh$nB, sh$nC, sh$ratio_B, sh$ratio_C, sh$decision),
sprintf("- Tanı VIF: AD_B %.1f, GBM_B %.1f", b4$models$AD_B$vif_diag, b4$models$GBM_B$vif_diag),
"", "## H2 (ikincil; Stage 3) — plan karar kuralı",
"| Analiz | Test 1 DU | Test 1 DD | Test 2 | Test 3 | Karar |", "|---|---|---|---|---|---|",
t3(s3$reproduction_original, "Orijinal (psödo-replike, yaşsız, Braak II)"),
t3(s3$plan_B6, "Plan: Braak VI vs 0, donör, yaş"),
t3(s3$deviation_B2, "Sapma: Braak II vs 0, donör, yaş"),
"", "## Stage 1",
sprintf("- Braak VI vs 0: donör %s, artık df %d, FDR<0,05 = %d, imza = %d gen",
        paste(s1$primary_B6$donors, collapse = " vs "), s1$primary_B6$resid_df,
        sum(s1$primary_B6$tt$FDR < 0.05), length(s1$signature_B6)),
sprintf("- Braak II vs 0: donör %s, artık df %d, FDR<0,05 = %d, imza = %d gen",
        paste(s1$primary_B2$donors, collapse = " vs "), s1$primary_B2$resid_df,
        sum(s1$primary_B2$tt$FDR < 0.05), length(s1$signature_B2)),
sprintf("- Psödo-replikasyon: donör ICC %.3f; FDR<0,05 orijinal %d, voom+dupCor %d, donör-toplam %d; logFC r = %.3f",
        st$pseudoreplication$icc_donor, st$pseudoreplication$n_sig["original"],
        st$pseudoreplication$n_sig["voom_dupcor"], st$pseudoreplication$n_sig["donor_sum"],
        st$pseudoreplication$logFC_r_orig_vs_donor),
"", "## Negatif kontrol (keşifsel)",
sprintf("- Kol 1 (gürültü): AD korunma %.3f -> %s", nc$retain_arm1, nc$arm1_call),
sprintf("- Kol 2 (rastgele gen PC): AD korunma %.3f, medyan kompozisyon R² %.3f -> %s",
        nc$retain_arm2, nc$median_compR2, nc$arm2_call),
"", "## GSE125583 (keşifsel)",
sprintf("- Yaş-grup r = %.3f; stres genleri yaş+cinsiyet+kompozisyon ayarlı: %d/%d yukarı, %d/%d FDR<0,05",
        st$validation_age_adjusted$r_group_age, sum(dd$logFC_comp[S] > 0), sum(S),
        sum(S & dd$logFC_comp > 0 & dd$FDR_comp < 0.05), sum(S)))
writeLines(facts, here::here("manuscript/facts_scenarioB.md"))
cat(facts, sep = "\n")

