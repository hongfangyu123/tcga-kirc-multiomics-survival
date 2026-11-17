# -----------------------------------------------------------------------------
# Part B3. Multivariable prognostic models: penalised Cox regression.
# Elastic net (alpha = 0.2) and LASSO (alpha = 1) on stage + age + the screened
# molecular features; evaluated on 100 random 50/50 train/test splits with the
# penalty re-tuned by cross-validation inside every training half.
# -----------------------------------------------------------------------------
source("R/00_setup.R")
suppressPackageStartupMessages(library(glmnet))

sf <- readRDS(file.path(OUT_DIR, "screened_features.rds"))   # from 05
s  <- sf$samples
s$stage3 <- factor(s$stage3, levels = c("Stage I/II", "Stage III", "Stage IV"))

# Design matrix: stage as dummies (its steps in risk are unequal: HR ~2 for III,
# ~8 for IV), age, and the molecular features with patients on the rows
clin_mm <- model.matrix(~ stage3 + age, data = s)[, -1]
colnames(clin_mm) <- c("StageIII", "StageIV", "Age")
Xmat <- cbind(clin_mm, t(sf$X))
rownames(Xmat) <- s$patient
Ysurv <- Surv(s$event.time, s$event)
cat("Design matrix:", nrow(Xmat), "patients x", ncol(Xmat), "variables\n")
saveRDS(list(Xmat = Xmat, Ysurv = Ysurv, samples = s, feat = sf$feat),
        file.path(OUT_DIR, "design.rds"))

# Elastic net on all patients, lambda tuned by 10-fold CV on the C-index -----------
set.seed(5241)
foldid  <- sample(rep(1:10, length.out = nrow(Xmat)))
cv_enet <- cv.glmnet(Xmat, Ysurv, family = "cox", alpha = 0.2,
                     type.measure = "C", foldid = foldid)
save_png("13_elastic_net_cv_curve.png", plot(cv_enet), width = 7, height = 5)
print(cv_enet)

cf <- as.matrix(coef(cv_enet, s = "lambda.min"))[, 1]
cf <- cf[cf != 0]
enet_tab <- data.frame(Variable = names(cf), coefficient = unname(cf))
ann <- sf$feat[, c("Variable", "VarName", "type", "gene", "chr", "pos", "relatedLocation")]
enet_tab <- merge(enet_tab, ann, by = "Variable", all.x = TRUE, sort = FALSE)
enet_tab$type[enet_tab$Variable %in% colnames(clin_mm)] <- "clinical"
enet_tab <- enet_tab[order(-abs(enet_tab$coefficient)), ]
cat("Non-zero coefficients:", nrow(enet_tab), "\n"); print(table(enet_tab$type))
write.csv(enet_tab, file.path(RES_DIR, "elastic_net_coefficients.csv"), row.names = FALSE)

# 100 random 50/50 splits (saved, so the random forest uses exactly the same) -------
n_rep <- 100
n <- nrow(Xmat)
split_file <- file.path(OUT_DIR, "splits.rds")
if (file.exists(split_file)) {
  splits <- readRDS(split_file)
} else {
  set.seed(2026)
  splits <- lapply(seq_len(n_rep), function(r) sort(sample(n, floor(n / 2))))
  saveRDS(splits, split_file)
}

res_file <- file.path(OUT_DIR, "penalised_cox_splits.rds")
if (file.exists(res_file)) {
  cv_res <- readRDS(res_file)
} else {
  set.seed(2026)
  cv_res <- data.frame(rep = 1:n_rep, C_lasso = NA, C_enet = NA, C_clin = NA,
                       nz_lasso = NA, nz_enet = NA)
  for (r in 1:n_rep) {
    tr  <- splits[[r]]
    te  <- setdiff(1:n, tr)
    fid <- sample(rep(1:10, length.out = length(tr)))
    for (m in c("lasso", "enet")) {
      a  <- if (m == "lasso") 1 else 0.2
      cv <- cv.glmnet(Xmat[tr, ], Ysurv[tr], family = "cox", alpha = a,
                      type.measure = "C", foldid = fid)      # tuned in training half only
      risk <- predict(cv, newx = Xmat[te, ], s = "lambda.min")
      cv_res[r, paste0("C_", m)]  <- Cindex(risk, Ysurv[te])
      cv_res[r, paste0("nz_", m)] <- cv$nzero[cv$index["min", 1]]
    }
    # Clinical-only Cox model on the same split: the fair benchmark
    fit_c  <- coxph(Surv(event.time, event) ~ age + stage3, data = s[tr, ])
    cv_res$C_clin[r] <- Cindex(predict(fit_c, newdata = s[te, ], type = "lp"), Ysurv[te])
    if (r %% 10 == 0) cat("split", r, "of", n_rep, "done\n")
  }
  saveRDS(cv_res, res_file)
}
print(summary(cv_res[, -1]))

save_png("14_lasso_vs_elastic_net.png",
         boxplot(cv_res$C_lasso, cv_res$C_enet, names = c("LASSO", "Elastic net (0.2)"),
                 col = c("orange", "skyblue"), ylab = "Held-out concordance index",
                 main = "100 random 50/50 splits"),
         width = 7, height = 5)
print(wilcox.test(cv_res$C_enet, cv_res$C_lasso, paired = TRUE))
