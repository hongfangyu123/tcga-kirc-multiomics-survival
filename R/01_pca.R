# -----------------------------------------------------------------------------
# Part A1. Unsupervised view: do the main axes of molecular variation relate to
# survival status?  PCA of each omics block, coloured by vital status.
# -----------------------------------------------------------------------------
source("R/00_setup.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

d <- load_kirc()
is_meth <- d$feat$type == "methylation"
is_rna  <- d$feat$type == "mRNA"
cat("Patients:", ncol(d$X), "| methylation:", sum(is_meth), "| mRNA:", sum(is_rna), "\n")

# Patients on rows; drop features with missing values or zero variance
prep_pca <- function(M) {
  M <- t(M)
  M <- M[, colSums(is.na(M)) == 0, drop = FALSE]
  M[, apply(M, 2, sd) > 0, drop = FALSE]
}
pca_m <- prcomp(prep_pca(d$X[is_meth, ]), scale. = TRUE, rank. = 4)
pca_r <- prcomp(prep_pca(d$X[is_rna, ]),  scale. = TRUE, rank. = 4)
ve_m <- 100 * pca_m$sdev^2 / sum(pca_m$sdev^2)
ve_r <- 100 * pca_r$sdev^2 / sum(pca_r$sdev^2)

ve_tab <- rbind(Methylation = ve_m[1:4], mRNA = ve_r[1:4])
colnames(ve_tab) <- paste0("PC", 1:4)
print(round(ve_tab, 2))
write.csv(round(ve_tab, 2), file.path(RES_DIR, "pca_variance_explained.csv"))

status <- factor(d$samples$vital_status)
pc_plot <- function(pca, ve, a, b, ttl) {
  df <- data.frame(x = pca$x[, a], y = pca$x[, b], Status = status)
  ggplot(df, aes(x, y, colour = Status)) +
    geom_point(size = 1.6, alpha = 0.7) +
    stat_ellipse(level = 0.95) +
    labs(title = ttl, x = sprintf("PC%d (%.1f%%)", a, ve[a]),
         y = sprintf("PC%d (%.1f%%)", b, ve[b])) +
    theme_bw()
}
p <- (pc_plot(pca_m, ve_m, 1, 2, "Methylation: PC1 vs PC2") |
      pc_plot(pca_m, ve_m, 2, 3, "Methylation: PC2 vs PC3")) /
     (pc_plot(pca_r, ve_r, 1, 2, "mRNA: PC1 vs PC2") |
      pc_plot(pca_r, ve_r, 1, 3, "mRNA: PC1 vs PC3")) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "PCA of TCGA-KIRC tumours by vital status")
save_png("01_pca_by_vital_status.png", p, width = 10, height = 8)
