# =============================================================================
# 07e_stage4_negative_control.R — EXPLORATORY negative control for H1 (see DEVIATIONS pre-declaration)
# =============================================================================
suppressPackageStartupMessages({ library(limma); library(BRETIGEA) })
ad  <- readRDS(here::here("results/intermediate/ad_bulk.rds"))
gbm <- readRDS(here::here("results/intermediate/gbm_bulk.rds"))
keep_g <- !is.na(gbm$meta$sex) & gbm$meta$sex != ""
gbm$meta <- droplevels(gbm$meta[keep_g, ]); gbm$expr <- gbm$expr[, gbm$meta$sample_id]
ad$meta$diagnosis <- relevel(factor(ad$meta$diagnosis), ref = "Control")
b4 <- readRDS(here::here("outputs/stage4_H1_bretigea.rds"))
bret_ad <- as.matrix(b4$bretigea$ad[colnames(ad$expr), ])

markers <- unique(BRETIGEA::markers_df_human_brain$markers)
pool_ad  <- setdiff(rownames(ad$expr),  markers)
pool_gbm <- setdiff(rownames(gbm$expr), markers)

n_deg <- function(expr, des) {
  f <- eBayes(lmFit(expr, des)); cf <- grep("^diagnosis", colnames(des), value = TRUE)[1]
  tt <- topTable(f, coef = cf, number = Inf, sort.by = "none"); tt$gene <- rownames(tt); tt
}
shared_n <- function(a, b) {
  g <- intersect(a$gene, b$gene); x <- a[match(g, a$gene), ]; y <- b[match(g, b$gene), ]
  sum(x$adj.P.Val < 0.05 & y$adj.P.Val < 0.05 & sign(x$logFC) == sign(y$logFC))
}
pc1 <- function(expr, genes) prcomp(t(expr[genes, ]), center = TRUE, scale. = TRUE)$x[, 1]

base_ad  <- model.matrix(~ diagnosis + sex + age + dataset, data = ad$meta)
base_gbm <- model.matrix(~ diagnosis + sex, data = gbm$meta)
N_IT <- 1000; set.seed(PARAMS$seed)
out <- vector("list", N_IT)
for (i in seq_len(N_IT)) {
  # Arm 1: noise
  z_ad  <- matrix(rnorm(nrow(base_ad) * 6),  ncol = 6, dimnames = list(NULL, paste0("z", 1:6)))
  z_gbm <- matrix(rnorm(nrow(base_gbm) * 6), ncol = 6, dimnames = list(NULL, paste0("z", 1:6)))
  a1 <- n_deg(ad$expr,  cbind(base_ad,  z_ad)); g1 <- n_deg(gbm$expr, cbind(base_gbm, z_gbm))
  # Arm 2: random-gene PCs
  r_ad  <- sapply(1:6, function(j) pc1(ad$expr,  sample(pool_ad, 50)));  colnames(r_ad)  <- paste0("r", 1:6)
  r_gbm <- sapply(1:6, function(j) pc1(gbm$expr, sample(pool_gbm, 50))); colnames(r_gbm) <- paste0("r", 1:6)
  a2 <- n_deg(ad$expr,  cbind(base_ad,  r_ad)); g2 <- n_deg(gbm$expr, cbind(base_gbm, r_gbm))
  r2 <- mean(apply(bret_ad, 2, function(y) summary(lm(y ~ r_ad))$r.squared))
  out[[i]] <- c(arm1_ad = sum(a1$adj.P.Val < 0.05), arm1_shared = shared_n(a1, g1),
                arm2_ad = sum(a2$adj.P.Val < 0.05), arm2_shared = shared_n(a2, g2), arm2_compR2 = r2)
  if (i %% 50 == 0) message(sprintf("iteration %d / %d", i, N_IT))
}
nc <- as.data.frame(do.call(rbind, out))

retain1 <- median(nc$arm1_ad) / 4449; retain2 <- median(nc$arm2_ad) / 4449; R2 <- median(nc$arm2_compR2)
arm1_call <- if (retain1 >= 0.80) "power loss does NOT explain collapse" else
             if (retain1 < 0.20) "collapse largely POWER LOSS" else "PARTIAL"
arm2_call <- if (R2 >= 0.25) "NON-DISCRIMINATING (random genes capture composition)" else
             if (retain2 < 0.20) "collapse NOT composition-specific" else
             if (retain2 >= 0.80) "COMPOSITION-SPECIFIC" else "PARTIAL"
stage4_negctrl <- list(iterations = nc, retain_arm1 = retain1, retain_arm2 = retain2,
                       median_compR2 = R2, arm1_call = arm1_call, arm2_call = arm2_call,
                       run_time = Sys.time())
saveRDS(stage4_negctrl, here::here("outputs/stage4_negative_control.rds"))
message("Saved outputs/stage4_negative_control.rds")

