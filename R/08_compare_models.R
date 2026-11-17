# -----------------------------------------------------------------------------
# Part B5. Head-to-head comparison on the same 100 splits, and how far the
# elastic net and the forest agree on which features matter.
# -----------------------------------------------------------------------------
source("R/00_setup.R")

cv_res <- readRDS(file.path(OUT_DIR, "penalised_cox_splits.rds"))
cv_res$C_rsf <- readRDS(file.path(OUT_DIR, "rsf_splits.rds"))
cindex_clin  <- readRDS(file.path(OUT_DIR, "cindex_clinical.rds"))
cmp <- cv_res[, c("C_lasso", "C_enet", "C_rsf", "C_clin")]

save_png("17_model_comparison.png", {
  boxplot(cmp$C_lasso, cmp$C_enet, cmp$C_rsf, names = c("LASSO", "Elastic net", "Random forest"),
          col = c("orange", "skyblue", "palegreen"), ylab = "Held-out concordance index",
          main = "Same 100 random 50/50 splits")
  abline(h = cindex_clin, col = "red", lty = 2, lwd = 2)
  abline(h = median(cmp$C_clin), col = "red", lty = 3, lwd = 2)
  legend("bottomright", bty = "n", lty = c(2, 3), col = "red", lwd = 2, cex = 0.8,
         legend = c("Clinical Cox, apparent C", "Clinical Cox, held-out median"))
}, width = 7, height = 5)

summ <- data.frame(model  = c("LASSO Cox", "Elastic net Cox", "Random survival forest",
                              "Clinical Cox (age + stage)"),
                   mean_C   = colMeans(cmp),
                   median_C = apply(cmp, 2, median),
                   sd_C     = apply(cmp, 2, sd))
summ[, -1] <- round(summ[, -1], 3)
print(summ, row.names = FALSE)
write.csv(summ, file.path(RES_DIR, "model_comparison.csv"), row.names = FALSE)

for (pr in list(c("C_enet", "C_lasso"), c("C_rsf", "C_enet"),
                c("C_enet", "C_clin"), c("C_lasso", "C_clin"), c("C_rsf", "C_clin"))) {
  a <- cmp[[pr[1]]]; b <- cmp[[pr[2]]]
  cat(pr[1], "vs", pr[2], ": mean diff", round(mean(a - b), 3),
      "| first better in", sum(a > b), "of 100 | paired Wilcoxon p =",
      signif(wilcox.test(a, b, paired = TRUE)$p.value, 2), "\n")
}

# Feature agreement between elastic net and forest ------------------------------------
dz   <- readRDS(file.path(OUT_DIR, "design.rds"))
rsf  <- readRDS(file.path(OUT_DIR, "rsf_full.rds"))
enet <- read.csv(file.path(RES_DIR, "elastic_net_coefficients.csv"))
rf_names <- make.names(colnames(dz$Xmat))          # names as used inside rfsrc
vimp_all <- sort(rsf$importance, decreasing = TRUE)
names(vimp_all) <- colnames(dz$Xmat)[match(names(vimp_all), rf_names)]
top_rsf <- names(vimp_all)[1:30]
cat("Forest top 30 also selected by the elastic net:",
    length(intersect(top_rsf, enet$Variable)), "\n")

# Correlated groups explain the disagreement: each top forest feature that the
# elastic net dropped is strongly correlated with one that it kept
mol <- dz$Xmat[, -(1:3)]
miss <- setdiff(intersect(top_rsf, colnames(mol)), enet$Variable)
kept <- intersect(enet$Variable, colnames(mol))
max_r <- sapply(miss, function(f) max(abs(cor(mol[, f], mol[, kept]))))
print(round(sort(max_r, decreasing = TRUE), 2))

# Methylation vs expression of the same gene: SLC22A2 and AQP1
for (g in c("SLC22A2", "AQP1")) {
  hit <- dz$feat$Variable[vapply(strsplit(dz$feat$gene, ";", fixed = TRUE),
                                 function(x) g %in% x, logical(1))]
  cat("\n", g, "- features passing the screen:", length(hit), "\n")
  if (length(hit)) {
    print(data.frame(feature = hit, type = dz$feat$type[match(hit, dz$feat$Variable)],
                     in_elastic_net = hit %in% enet$Variable,
                     rsf_rank = match(hit, names(vimp_all))), row.names = FALSE)
    if (length(hit) > 1) print(round(cor(mol[, hit]), 2))
  }
}
