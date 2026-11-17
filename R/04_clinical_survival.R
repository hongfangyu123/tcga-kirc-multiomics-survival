# -----------------------------------------------------------------------------
# Part B1. Clinical benchmark: what do stage, age and gender already tell us?
# Kaplan-Meier curves, log-rank / Peto-Peto tests and a clinical Cox model.
# Patients without AJCC stage are excluded from here on.
# -----------------------------------------------------------------------------
source("R/00_setup.R")
suppressPackageStartupMessages({ library(survminer); library(ggplot2) })

d <- staged_cohort(load_kirc())
s <- d$samples
cat("Patients with stage:", nrow(s), "| deaths:", sum(s$event), "\n")
y <- Surv(s$event.time, s$event)

km_plot <- function(fit, file, title, labs = NULL, ...) {
  args <- list(fit, data = s, risk.table = TRUE, break.time.by = 730,
               xlab = "Days since diagnosis", title = title, ggtheme = theme_bw(), ...)
  if (!is.null(labs)) args <- c(args, list(legend.labs = labs, legend.title = ""))
  save_png(file, do.call(ggsurvplot, args), width = 6.5, height = 6.5)
}

# Overall survival
km_all <- survfit(y ~ 1)
print(km_all)
print(summary(km_all, times = c(365, 1095, 1825)))
km_plot(km_all, "06_km_overall.png", sprintf("Overall survival (n = %d)", nrow(s)),
        conf.int = TRUE, surv.median.line = "hv", legend = "none")

# Stage: four groups, then I and II merged (their curves cross; log-rank p = 0.35)
s$stage4 <- factor(s$stage, 1:4, c("Stage I", "Stage II", "Stage III", "Stage IV"))
km_plot(survfit(Surv(event.time, event) ~ stage4, data = s), "07_km_by_stage_4groups.png",
        "Survival by stage (4 groups)", levels(s$stage4), pval = TRUE)
for (pr in list(1:2, 2:3, 3:4)) {
  sd <- survdiff(Surv(event.time, event) ~ stage, data = s, subset = stage %in% pr)
  cat("Stage", pr[1], "vs", pr[2], ": log-rank p =", signif(1 - pchisq(sd$chisq, 1), 2), "\n")
}
km_plot(survfit(Surv(event.time, event) ~ stage3, data = s), "08_km_by_stage_3groups.png",
        "Survival by stage (I and II merged)", levels(s$stage3), pval = TRUE)

# Log-rank (rho = 0, equal weights) vs Peto-Peto (rho = 1, early events weighted more)
print(survdiff(y ~ s$stage3, rho = 0))
print(survdiff(y ~ s$stage3, rho = 1))

# Gender and age group
km_plot(survfit(Surv(event.time, event) ~ gender, data = s), "09_km_by_gender.png",
        "Survival by gender", c("Female", "Male"), pval = TRUE)
s$age_group <- factor(ifelse(s$age >= 60, ">= 60", "< 60"), levels = c("< 60", ">= 60"))
km_plot(survfit(Surv(event.time, event) ~ age_group, data = s), "10_km_by_age.png",
        "Survival by age group", c("Age < 60", "Age >= 60"), pval = TRUE)

# Clinical Cox model: the benchmark for every molecular model
cox_clin <- coxph(y ~ age + stage3 + gender, data = s)
print(summary(cox_clin))
cox_clin2 <- coxph(y ~ age + stage3, data = s)
print(anova(cox_clin2, cox_clin))          # gender adds nothing (p = 0.97)

hr <- summary(cox_clin)$conf.int
clin_tab <- data.frame(term = rownames(hr), HR = hr[, 1], lower95 = hr[, 3], upper95 = hr[, 4],
                       p = summary(cox_clin)$coefficients[, 5])
write.csv(clin_tab, file.path(RES_DIR, "clinical_cox_model.csv"), row.names = FALSE)
cindex_clin <- unname(cox_clin$concordance["concordance"])
cat("Clinical model C-index (apparent):", round(cindex_clin, 3), "\n")
saveRDS(cindex_clin, file.path(OUT_DIR, "cindex_clinical.rds"))
