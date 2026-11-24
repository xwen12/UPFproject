# Total UPF and 8 UPF subclasses vs the 42 UPF-related lipids (same mixed model as 2.1).
# Input: analysis datasets (1.3), 2.1 results. Output: result/UPFsubcat_lipids42_adj.csv.

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
raw_dir  <- "path/to/SHS_raw_data"

library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)

rawData   <- raw_dir
readyData <- file.path(proj_dir, "data")
result    <- file.path(proj_dir, "result")
tabsFigs  <- file.path(proj_dir, "result", "tabsFigs")

setwd(readyData)
load("lipidomics_p4.Rda")
load("lipidomics_p5.Rda")
load("dat_AimOne_p4.Rda")
load("dat_AimOne_p5.Rda")

dat_AimOne_p4 <- dat_AimOne_p4 %>%
  select(-chd, -cvd, -TIME_chd, -TIME_cvd, -starts_with("incident"))

temp_1 <- as.data.frame(lipidomics_p4) %>%
  inner_join(dat_AimOne_p4, by = c("ID" = "IDNO"))

temp_2 <- as.data.frame(lipidomics_p5) %>%
  inner_join(dat_AimOne_p5, by = c("ID" = "IDNO"))

temp_1[, "Phase"] <- "Phase4"
temp_2[, "Phase"] <- "Phase5"
names(temp_1) <- sub("^S4", "", names(temp_1))
names(temp_2) <- sub("^S5", "", names(temp_2))
colnames(temp_1)[colnames(temp_1) == "PDM3"]        <- "PA"
colnames(temp_2)[colnames(temp_2) == "ped_avg3day"] <- "PA"
temp_2$BOOKNUM <- NULL

temp_3 <- rbind(temp_1, temp_2)
temp_3$Phase <- as.factor(temp_3$Phase)

dat <- temp_3 %>%
  mutate(across(c(ends_with("_P"), AGE, EDU, BMI, PA, DT_KCAL),
                ~ as.numeric(scale(.x, center = TRUE, scale = TRUE))))

nova_vars <- c("nova4_kcal_P", paste0("nova4_", 1:8, "_kcal_P"))

dt        <- dat
first_ind <- which(colnames(dt) == "AC(10:0)")
last_ind  <- which(colnames(dt) == "Unknown_1230")
for (i in first_ind:last_ind) {
  dt[, i] <- as.numeric(scale(dt[, i], center = TRUE, scale = TRUE))
}

dt[, "FID"] <- as.factor(dt[, "FID"])
dt[, "ID"]  <- as.factor(dt[, "ID"])

setwd(result)
UPF_lipids_all <- read.table("UPF_lipids_all.csv", sep = ",", header = TRUE,
                             check.names = FALSE)
selected_lipids_names <- UPF_lipids_all$lipids_name

covars <- paste("AGE + SEX + CENTER + EDU + BMI + SMOKE + PA + ETOH + LIPIDLOW +",
                "DT_KCAL + Phase + (1|FID) + (1|ID)")

results_list <- list()
log_list     <- list()

for (i in seq_along(selected_lipids_names)) {
  lipid_name   <- selected_lipids_names[i]
  dt$lipid_tmp <- dt[[lipid_name]]

  for (nova in nova_vars) {
    fml  <- as.formula(paste0("lipid_tmp ~ ", nova, " + ", covars))
    warn <- character(0)

    fit <- tryCatch(
      withCallingHandlers(
        suppressMessages(lmer(fml, data = dt)),
        warning = function(w) {
          warn <<- c(warn, conditionMessage(w))
          invokeRestart("muffleWarning")
        }),
      error = function(e) { warn <<- c(warn, paste("ERROR:", conditionMessage(e))); NULL }
    )

    if (is.null(fit)) {
      results_list[[length(results_list) + 1]] <- data.frame(
        lipids_name = lipid_name, nova = nova,
        beta = NA_real_, se = NA_real_, lower_CI = NA_real_, upper_CI = NA_real_,
        pvals = NA_real_, n_obs = NA_integer_, n_id = NA_integer_)
      log_list[[length(log_list) + 1]] <- data.frame(
        lipids_name = lipid_name, nova = nova, failed = TRUE, singular = NA,
        warnings = paste(warn, collapse = " | "))
      next
    }

    a <- coef(summary(fit))
    b <- a[nova, "Estimate"]
    s <- a[nova, "Std. Error"]

    results_list[[length(results_list) + 1]] <- data.frame(
      lipids_name = lipid_name, nova = nova,
      beta = b, se = s, lower_CI = b - 1.96 * s, upper_CI = b + 1.96 * s,
      pvals = a[nova, "Pr(>|t|)"],
      n_obs = nobs(fit), n_id = as.integer(ngrps(fit)["ID"]))

    log_list[[length(log_list) + 1]] <- data.frame(
      lipids_name = lipid_name, nova = nova, failed = FALSE,
      singular = isSingular(fit, tol = 1e-4),
      warnings = paste(warn, collapse = " | "))
  }
}

res     <- bind_rows(results_list)
fit_log <- bind_rows(log_list)

res_sub <- res %>%
  filter(nova != "nova4_kcal_P") %>%
  mutate(qvals = p.adjust(pvals, method = "BH"))

res_tot <- res %>%
  filter(nova == "nova4_kcal_P") %>%
  left_join(UPF_lipids_all %>% select(lipids_name, qvals), by = "lipids_name")

res <- bind_rows(res_tot, res_sub) %>%
  mutate(fdr_group = ifelse(nova == "nova4_kcal_P", "total_1542", "subcat_336"))

res <- res %>%
  mutate(sig = case_when(qvals < 0.05 ~ "**",
                         pvals < 0.05 ~ "*",
                         TRUE         ~ ""))

nova_labels <- c(
  "nova4_kcal_P"   = "Total UPF",
  "nova4_1_kcal_P" = "Ultra-processed breads & breakfast foods",
  "nova4_2_kcal_P" = "Sauces, spreads & condiments",
  "nova4_3_kcal_P" = "Packaged sweet snacks & desserts",
  "nova4_4_kcal_P" = "Packaged savory snacks",
  "nova4_5_kcal_P" = "Sugar & artificially sweetened beverages",
  "nova4_6_kcal_P" = "Meat/poultry/seafood ready-to-eat products",
  "nova4_7_kcal_P" = "Ready-to-eat/heat mixed dishes",
  "nova4_8_kcal_P" = "Dairy-based desserts"
)
res$nova_label <- unname(nova_labels[res$nova])

setwd(result)
write.table(res,     "UPFsubcat_lipids42_adj.csv",        sep = ",", row.names = F, col.names = T)
write.table(fit_log, "UPFsubcat_lipids42_adj_fitlog.csv", sep = ",", row.names = F, col.names = T)
