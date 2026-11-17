# -----------------------------------------------------------------------------
# Shared settings and helper functions used by every analysis script.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

# Input folder with the three processed TCGA-KIRC tables (see data/README.md).
# Override with:  Sys.setenv(KIRC_DATA_DIR = "/path/to/folder")
DATA_DIR <- Sys.getenv("KIRC_DATA_DIR", "data")
OUT_DIR  <- "output"    # large intermediate objects (git-ignored)
FIG_DIR  <- "figures"   # figures shown in the README
RES_DIR  <- "results"   # small result tables kept in the repository
for (d in c(OUT_DIR, FIG_DIR, RES_DIR)) dir.create(d, showWarnings = FALSE)

FDR <- 0.01   # false discovery rate used throughout

# Save a base-R or ggplot figure as PNG -----------------------------------------
save_png <- function(file, expr, width = 8, height = 6, res = 150) {
  png(file.path(FIG_DIR, file), width = width, height = height,
      units = "in", res = res)
  on.exit(dev.off())
  out <- expr
  # ggplot-type objects are drawn only when printed; base graphics and
  # ComplexHeatmap::draw() have already drawn by this point
  if (inherits(out, c("gg", "ggplot", "patchwork", "ggsurvplot")))
    print(out)
  invisible(out)
}

# Load and prepare the integrated data -------------------------------------------
# Returns a list with
#   X       numeric matrix, molecular features (rows) x patients (columns)
#   feat    feature annotation, row-aligned with X
#   samples patient table, row-aligned with the columns of X
load_kirc <- function(data_dir = DATA_DIR) {
  f <- function(x) file.path(data_dir, x)
  stopifnot(file.exists(f("integrated_data.txt")))

  # fread is much faster than read.delim for the ~750 MB matrix;
  # check.names = FALSE keeps TCGA barcodes (which contain "-") intact
  mdata   <- fread(f("integrated_data.txt"), data.table = FALSE, check.names = FALSE)
  feat    <- fread(f("variable_info.txt"),   data.table = FALSE, check.names = FALSE)
  samples <- fread(f("sample_info.txt"),     data.table = FALSE, check.names = FALSE)

  stopifnot(identical(feat$Variable, mdata$Variable))   # rows aligned

  # Keep DNA methylation and mRNA only
  keep <- mdata$Type %in% c("methylation", "mRNA")
  mdata <- mdata[keep, ]
  feat  <- feat[keep, ]

  # Annotation clean-up: relatedGene repeats symbols ("A;A;B"); keep the raw
  # string as well, because its positions pair one-to-one with relatedLocation
  feat$relatedGene[is.na(feat$relatedGene)] <- ""
  feat$relatedLocation[is.na(feat$relatedLocation)] <- ""
  feat$relatedGene_raw <- feat$relatedGene
  feat$relatedGene <- vapply(strsplit(feat$relatedGene, ";", fixed = TRUE),
                             function(g) paste(unique(g), collapse = ";"),
                             character(1))
  # One gene label per feature: related gene(s) for CpG probes, symbol for mRNA
  feat$gene <- ifelse(feat$type == "mRNA", feat$VarName, feat$relatedGene)

  # Survival outcome: event = 1 if deceased, 0 if censored
  samples$event      <- as.integer(samples$vital_status == "Dead")
  samples$event.time <- as.numeric(samples$days_to_last_follow_up)
  samples <- samples[!is.na(samples$event.time) & samples$event.time > 0, ]
  samples$age <- samples$age_at_index

  # AJCC stage -> ordinal 1-4 via a lookup table, so that a missing or
  # unexpected label stays NA instead of silently falling into a category
  stage_map <- c("Stage I" = 1, "Stage IA" = 1, "Stage IB" = 1,
                 "Stage II" = 2, "Stage III" = 3, "Stage IV" = 4)
  samples$stage <- unname(stage_map[samples$ajcc_pathologic_stage])
  # Stage I and II have indistinguishable survival (see 04_clinical_survival.R)
  samples$stage3 <- factor(c("Stage I/II", "Stage I/II", "Stage III", "Stage IV")[samples$stage],
                           levels = c("Stage I/II", "Stage III", "Stage IV"))

  # Align the matrix columns with the patient table
  col_idx <- match(samples$patient, colnames(mdata))
  stopifnot(!anyNA(col_idx))
  X <- as.matrix(mdata[, col_idx])
  rownames(X) <- mdata$Variable
  rm(mdata); invisible(gc())

  list(X = X, feat = feat, samples = samples)
}

# Patients with recorded stage, ordered by stage then follow-up time
staged_cohort <- function(d) {
  keep <- which(!is.na(d$samples$stage))
  keep <- keep[order(d$samples$stage[keep], d$samples$event.time[keep])]
  d$samples <- d$samples[keep, ]
  rownames(d$samples) <- NULL
  d$X <- d$X[, keep, drop = FALSE]
  d
}
