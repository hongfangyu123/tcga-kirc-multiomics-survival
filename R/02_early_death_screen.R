# -----------------------------------------------------------------------------
# Part A2. Supervised case-control screen.
# Cases:    died within 2 years of diagnosis.
# Controls: alive and still under follow-up beyond 2 years.
# Patients censored before 2 years (outcome unknown) and patients who died after
# 2 years are excluded, so the contrast is clean.
# Every CpG probe and every gene is tested with a Welch two-sample t-test;
# FDR is controlled per platform with Storey's q-value.
# -----------------------------------------------------------------------------
source("R/00_setup.R")
suppressPackageStartupMessages({ library(qvalue); library(ggplot2); library(patchwork) })

d <- load_kirc()
s <- d$samples
case_idx <- which(s$event == 1 & s$event.time <= 730)
ctrl_idx <- which(s$event == 0 & s$event.time > 730)
cat("Cases:", length(case_idx), "| controls:", length(ctrl_idx),
    "| excluded:", nrow(s) - length(case_idx) - length(ctrl_idx), "\n")

# Welch t-test for every row of a block, with a progress counter
t_screen <- function(M, label, every) {
  n <- nrow(M)
  pval <- rep(NA_real_, n)
  diff <- rep(NA_real_, n)
  for (j in seq_len(n)) {
    x1 <- M[j, case_idx]; x1 <- x1[!is.na(x1)]
    x0 <- M[j, ctrl_idx]; x0 <- x0[!is.na(x0)]
    if (length(x1) >= 3 && length(x0) >= 3) {
      diff[j] <- mean(x1) - mean(x0)                 # case minus control
      tt <- tryCatch(t.test(x1, x0), error = function(e) NULL)
      if (!is.null(tt)) pval[j] <- tt$p.value
    }
    if (j %% every == 0) cat(label, j, "/", n, "\n")
  }
  data.frame(Variable = rownames(M), diff = diff, p = pval)
}

cache <- file.path(OUT_DIR, "early_death_ttests.rds")
if (file.exists(cache)) {
  tt <- readRDS(cache)
} else {
  is_meth <- d$feat$type == "methylation"
  tt <- rbind(cbind(t_screen(d$X[is_meth, ], "methylation", 25000), type = "methylation"),
              cbind(t_screen(d$X[!is_meth, ], "mRNA", 2000),        type = "mRNA"))
  tt <- tt[match(rownames(d$X), tt$Variable), ]
  saveRDS(tt, cache)
}

# p-value histograms ------------------------------------------------------------
br <- seq(0, 1, length.out = 101)
save_png("02_pvalue_histograms_early_death.png", {
  par(mfrow = c(1, 3))
  hist(tt$p[tt$type == "methylation"], breaks = br, main = "Methylation", xlab = "p-value")
  hist(tt$p[tt$type == "mRNA"],        breaks = br, main = "mRNA",        xlab = "p-value")
  hist(tt$p,                           breaks = br, main = "Pooled",      xlab = "p-value")
}, width = 11, height = 3.8)

# q-values: per platform vs pooled -------------------------------------------------
ok <- !is.na(tt$p)
tt$q_sep <- NA_real_
for (ty in c("methylation", "mRNA")) {
  k <- ok & tt$type == ty
  tt$q_sep[k] <- qvalue(tt$p[k])$qvalues
}
tt$q_pool <- NA_real_
tt$q_pool[ok] <- qvalue(tt$p[ok])$qvalues

sig_tab <- data.frame(
  strategy    = c("Per platform", "Pooled"),
  methylation = c(sum(tt$q_sep[tt$type == "methylation"]  < FDR, na.rm = TRUE),
                  sum(tt$q_pool[tt$type == "methylation"] < FDR, na.rm = TRUE)),
  mRNA        = c(sum(tt$q_sep[tt$type == "mRNA"]  < FDR, na.rm = TRUE),
                  sum(tt$q_pool[tt$type == "mRNA"] < FDR, na.rm = TRUE)))
print(sig_tab)
write.csv(sig_tab, file.path(RES_DIR, "early_death_fdr_by_strategy.csv"), row.names = FALSE)

save_png("03_qvalues_separate_vs_pooled.png", {
  par(mfrow = c(1, 2))
  for (ty in c("methylation", "mRNA")) {
    k <- tt$type == ty
    plot(tt$q_pool[k], tt$q_sep[k], pch = 20, cex = 0.3, main = ty,
         xlab = "q-value (pooled)", ylab = "q-value (per platform)")
    abline(0, 1, col = "red")
  }
}, width = 10, height = 4.5)

# Volcano plots (per-platform q-values) ---------------------------------------------
# Effect thresholds: 0.2 in beta-value for methylation, log2(1.5) for mRNA
thr <- c(methylation = 0.2, mRNA = log2(1.5))
tt$status <- "Not significant"
sig <- !is.na(tt$q_sep) & tt$q_sep < FDR
tt$status[sig & tt$diff >  thr[tt$type]] <- "Up"
tt$status[sig & tt$diff < -thr[tt$type]] <- "Down"
tt$status <- factor(tt$status, levels = c("Down", "Not significant", "Up"))
print(table(tt$type, tt$status))

volcano <- function(ty) {
  df <- tt[tt$type == ty & ok, ]
  # p-value matching q = 0.01 depends on the data, so it is read off the q-values
  pcut <- max(df$p[df$q_sep < FDR])
  ggplot(df, aes(diff, -log10(p), colour = status)) +
    geom_point(size = 0.4, alpha = 0.6) +
    scale_colour_manual(values = c(Down = "blue", "Not significant" = "grey70", Up = "red"),
                        drop = FALSE) +
    geom_vline(xintercept = c(-1, 1) * thr[[ty]], linetype = "dashed") +
    geom_hline(yintercept = -log10(pcut), linetype = "dashed") +
    labs(title = ty, x = "Mean difference (early death - long survivor)",
         y = "-log10(p-value)", colour = NULL) +
    theme_bw()
}
save_png("04_volcano_early_death.png",
         (volcano("methylation") | volcano("mRNA")) + plot_layout(guides = "collect"),
         width = 11, height = 5)

# Significant features, annotated -------------------------------------------------------
sig_out <- cbind(tt[tt$status != "Not significant", ],
                 d$feat[tt$status != "Not significant", c("VarName", "gene", "chr", "pos", "relatedLocation")])
sig_out <- sig_out[order(sig_out$p), ]
write.csv(sig_out, file.path(RES_DIR, "early_death_significant_features.csv"), row.names = FALSE)
saveRDS(tt, file.path(OUT_DIR, "early_death_results.rds"))
