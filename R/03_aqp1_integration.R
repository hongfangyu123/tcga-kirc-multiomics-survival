# -----------------------------------------------------------------------------
# Part A3. Do methylation and expression point to the same genes?
# Genes implicated on both platforms, then a close look at AQP1: does promoter
# methylation move opposite to its expression, as gene silencing predicts?
# -----------------------------------------------------------------------------
source("R/00_setup.R")
suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize) })

d  <- load_kirc()
tt <- readRDS(file.path(OUT_DIR, "early_death_results.rds"))   # from 02_early_death_screen.R
stopifnot(identical(tt$Variable, rownames(d$X)))
f <- d$feat

sig    <- tt$status != "Not significant"
is_met <- f$type == "methylation"
genes_m <- setdiff(unique(unlist(strsplit(f$relatedGene[sig & is_met], ";", fixed = TRUE))), "")
genes_r <- unique(f$VarName[sig & !is_met])
overlap <- sort(intersect(genes_m, genes_r))
cat("Genes from methylation:", length(genes_m), "| from mRNA:", length(genes_r),
    "| on both platforms:", length(overlap), "\n")
print(overlap)
writeLines(overlap, file.path(RES_DIR, "genes_significant_on_both_platforms.txt"))

# Target gene ----------------------------------------------------------------------
target <- "AQP1"
# grep("AQP1") would also match AQP10, AQP11, AQP12B ... so match whole symbols only
gene_lists <- strsplit(f$relatedGene, ";", fixed = TRUE)
meth_hits  <- which(is_met & vapply(gene_lists, function(g) target %in% g, logical(1)))
rna_hit    <- which(!is_met & f$VarName == target)
cat(target, ": 1 mRNA feature and", length(meth_hits), "CpG probes\n")

# Location of each probe relative to the target gene (raw gene / location
# strings are paired position by position)
loc <- vapply(meth_hits, function(k) {
  g <- strsplit(f$relatedGene_raw[k], ";", fixed = TRUE)[[1]]
  l <- strsplit(f$relatedLocation[k], ";", fixed = TRUE)[[1]]
  if (length(g) == length(l)) paste(unique(l[g == target]), collapse = ";") else f$relatedLocation[k]
}, character(1))

# Case / control patients as in 02_early_death_screen.R, cases first
s <- d$samples
case_idx <- which(s$event == 1 & s$event.time <= 730)
ctrl_idx <- which(s$event == 0 & s$event.time > 730)
ord <- c(case_idx, ctrl_idx)

M <- d$X[c(rna_hit, meth_hits), ord, drop = FALSE]
rownames(M) <- c(paste(target, "mRNA"), paste0(f$VarName[meth_hits], " (", loc, ")"))
Z <- t(scale(t(M)))                                  # z-score each feature
Z <- Z[rowSums(is.na(Z)) == 0, , drop = FALSE]

set.seed(5241)
grp <- factor(c(rep("Dead within 2y", length(case_idx)), rep("Alive beyond 2y", length(ctrl_idx))))
ha  <- HeatmapAnnotation(Outcome = grp,
                         col = list(Outcome = c("Dead within 2y" = "black", "Alive beyond 2y" = "gold")))
hm <- Heatmap(Z, name = "z-score", top_annotation = ha,
              col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
              clustering_distance_rows = "pearson", clustering_method_rows = "average",
              clustering_distance_columns = "pearson", clustering_method_columns = "average",
              show_column_names = FALSE, row_names_gp = gpar(fontsize = 8),
              column_title = paste(target, ": expression and methylation"))
save_png("05_aqp1_methylation_expression_heatmap.png", draw(hm), width = 10, height = 6)
