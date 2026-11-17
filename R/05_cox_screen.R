# -----------------------------------------------------------------------------
# Part B2. Genome-wide Cox screen.
# For each of ~346,000 features fit two Cox models for overall survival:
#   adjusted:    age + stage + feature
#   unadjusted:  feature only
# and keep the feature's coefficient, Wald p-value and the model C-index.
# Runtime: ~15 minutes on a laptop; the result is cached in output/.
# -----------------------------------------------------------------------------
source("R/00_setup.R")
suppressPackageStartupMessages(library(qvalue))

d <- staged_cohort(load_kirc())
s <- d$samples
y   <- Surv(s$event.time, s$event)
age <- s$age
st  <- s$stage3
cindex_clin <- readRDS(file.path(OUT_DIR, "cindex_clinical.rds"))   # from 04

cache <- file.path(OUT_DIR, "cox_screen.rds")
if (file.exists(cache)) {
  screen <- readRDS(cache)
} else {
  n_feat <- nrow(d$X)
  coef_adj <- p_adj <- c_adj <- coef_un <- p_un <- c_un <- rep(NA_real_, n_feat)
  for (i in seq_len(n_feat)) {
    x <- d$X[i, ]
    fa <- coxph(y ~ age + st + x)
    fu <- coxph(y ~ x)
    coef_adj[i] <- coef(fa)["x"]
    p_adj[i]    <- summary(fa)$coefficients["x", "Pr(>|z|)"]
    c_adj[i]    <- fa$concordance["concordance"]
    coef_un[i]  <- coef(fu)["x"]
    p_un[i]     <- summary(fu)$coefficients["x", "Pr(>|z|)"]
    c_un[i]     <- fu$concordance["concordance"]
    if (i %% 25000 == 0) cat(i, "of", n_feat, "features done\n")
  }
  screen <- data.frame(Variable = rownames(d$X), Type = d$feat$type,
                       coef_adj, p_adj, c_adj, coef_un, p_un, c_un)
  saveRDS(screen, cache)
}

# Note: low-variance CpG probes give huge coefficients (one unit = the whole
# 0-1 beta range), but the Wald z = beta/SE is scale-invariant, so p-values
# remain valid. Coefficients are therefore not used to rank features.

save_png("11_cox_screen_pvalue_histograms.png", {
  par(mfrow = c(1, 2))
  hist(screen$p_adj, breaks = 500, col = "steelblue", border = NA,
       main = "Adjusted (age + stage + feature)", xlab = "p-value")
  hist(screen$p_un, breaks = 500, col = "firebrick", border = NA,
       main = "Unadjusted (feature only)", xlab = "p-value")
}, width = 10, height = 5)

screen$q_adj <- qvalue(screen$p_adj)$qvalues
screen$q_un  <- qvalue(screen$p_un)$qvalues
cat("q <= 0.01: adjusted", sum(screen$q_adj <= FDR), "| unadjusted", sum(screen$q_un <= FDR), "\n")
print(table(Type = screen$Type, adjusted_q01 = screen$q_adj <= FDR))

save_png("12_qvalue_vs_cindex.png", {
  par(mfrow = c(1, 2))
  plot(screen$c_adj, -log10(screen$q_adj), pch = ".", col = rgb(0, 0, 0.6, 0.3),
       xlab = "C-index (age + stage + feature)", ylab = "-log10(q-value)", main = "Adjusted")
  abline(h = 2, col = "red", lty = 2); abline(v = cindex_clin, col = "darkgreen", lty = 2)
  plot(screen$c_un, -log10(screen$q_un), pch = ".", col = rgb(0.6, 0, 0, 0.3),
       xlab = "C-index (feature only)", ylab = "-log10(q-value)", main = "Unadjusted")
  abline(h = 2, col = "red", lty = 2); abline(v = 0.65, col = "darkgreen", lty = 2)
}, width = 10, height = 5)

# Full table (large, kept out of git) and the top 1,000 features (in git)
res <- cbind(screen[, 1:2], d$feat[, c("VarName", "gene")], screen[, -(1:2)])
res <- res[order(res$p_adj), ]
fwrite(res, file.path(OUT_DIR, "cox_screen_all_features.csv.gz"))
fwrite(head(res, 1000), file.path(RES_DIR, "cox_screen_top1000.csv"))

# Feature selection for the multivariable models:
# unadjusted C >= 0.65 (stand-alone discrimination) AND adjusted q <= 0.01
# (information beyond age and stage)
keep <- which(screen$c_un >= 0.65 & screen$q_adj <= FDR)
cat("Features passing the screen:", length(keep), "\n")
print(table(screen$Type[keep]))
saveRDS(list(X = d$X[keep, ], feat = d$feat[keep, ], samples = s),
        file.path(OUT_DIR, "screened_features.rds"))
