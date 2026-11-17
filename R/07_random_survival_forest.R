# -----------------------------------------------------------------------------
# Part B4. Random survival forest on the same design matrix: Brier score,
# subsampled variable importance with confidence intervals, and the same
# 100 train/test splits as the penalised Cox models.
# Memory-hungry: the 100-subsample importance step needs ~8 GB RAM.
# -----------------------------------------------------------------------------
source("R/00_setup.R")
suppressPackageStartupMessages({ library(randomForestSRC); library(glmnet) })

dz <- readRDS(file.path(OUT_DIR, "design.rds"))   # from 06
Xmat <- dz$Xmat; Ysurv <- dz$Ysurv; s <- dz$samples; feat <- dz$feat

# rfsrc may alter column names, so keep a map back to the original names
rf_df <- data.frame(event.time = s$event.time, event = s$event, Xmat)
name_map <- setNames(colnames(Xmat), colnames(rf_df)[-(1:2)])

rsf_file <- file.path(OUT_DIR, "rsf_full.rds")
if (file.exists(rsf_file)) {
  rsf <- readRDS(rsf_file)
} else {
  set.seed(5241)
  rsf <- rfsrc(Surv(event.time, event) ~ ., data = rf_df, ntree = 500, importance = TRUE)
  saveRDS(rsf, rsf_file)
}
print(rsf)

# Out-of-bag Brier score (censoring handled with Kaplan-Meier IPCW weights)
bs <- get.brier.survival(rsf, cens.model = "km")
cat("Integrated Brier score:", round(bs$crps.std, 3), "\n")
save_png("15_rsf_brier_score.png",
         plot(bs$brier.score, type = "s", col = 2, lwd = 2, xlab = "Days",
              ylab = "OOB Brier score", main = "Brier score over time"),
         width = 7, height = 5)

# Variable importance with delete-d jackknife CIs from 100 subsamples
sub_file <- file.path(OUT_DIR, "rsf_subsample.rds")
if (file.exists(sub_file)) {
  rsf_sub <- readRDS(sub_file)
} else {
  set.seed(5241)
  rsf_sub <- subsample(rsf, B = 100, verbose = FALSE)
  saveRDS(rsf_sub, sub_file)
}
vimp <- extract.subsample(rsf_sub)$var.jk.sel.Z
vimp$feature <- name_map[rownames(vimp)]
vimp <- vimp[order(-vimp$mean), ]
save_png("16_rsf_variable_importance.png",
         plot.subsample(rsf_sub, alpha = 0.05, xvar.names = rownames(vimp)[1:60],
                        xlab = "Variable importance (x 100)", cex = 0.6),
         width = 7, height = 14)

sig_vimp <- vimp[vimp$lower > 0, ]
sig_vimp$type <- ifelse(sig_vimp$feature %in% c("StageIII", "StageIV", "Age"), "clinical",
                        feat$type[match(sig_vimp$feature, feat$Variable)])
sig_vimp$gene <- feat$gene[match(sig_vimp$feature, feat$Variable)]
cat("Variables whose 95% CI lies above zero:", nrow(sig_vimp), "\n")
write.csv(sig_vimp[, c("feature", "type", "gene", "mean", "lower", "upper", "pvalue")],
          file.path(RES_DIR, "rsf_important_variables.csv"), row.names = FALSE)

# Same 100 splits as 06_penalised_cox.R
splits <- readRDS(file.path(OUT_DIR, "splits.rds"))
cache  <- file.path(OUT_DIR, "rsf_splits.rds")
if (file.exists(cache)) {
  C_rsf <- readRDS(cache)
} else {
  set.seed(2026)
  n <- nrow(Xmat)
  C_rsf <- rep(NA_real_, length(splits))
  for (r in seq_along(splits)) {
    tr <- splits[[r]]; te <- setdiff(1:n, tr)
    fit  <- rfsrc(Surv(event.time, event) ~ ., data = rf_df[tr, ], ntree = 500)
    pred <- predict(fit, newdata = rf_df[te, ])$predicted      # ensemble mortality
    C_rsf[r] <- Cindex(pred, Ysurv[te])
    if (r %% 10 == 0) cat("split", r, "of", length(splits), "done\n")
  }
  saveRDS(C_rsf, cache)
}
cat("RSF held-out C: mean", round(mean(C_rsf), 3), "| SD", round(sd(C_rsf), 3), "\n")
