# Sensitivity analysis: discovery in Arizona + South Dakota (P < 0.05), replication in Oklahoma
# (same sign, P < 0.05), same mixed model as 2.1; fixed-effect meta-analysis of the two sets.
# Output: result/sensCenter_AZSD_OK_*.csv.

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
raw_dir  <- "path/to/SHS_raw_data"
library(dplyr)
library(tidyr)
library(stringr)
library(lme4)
library(lmerTest)

rawData <- raw_dir
readyData <- file.path(proj_dir, "data")
result <-file.path(proj_dir, "result")

p_disc_cut  <- 0.05
p_repl_cut  <- 0.05
disc_sites  <- c("AZ", "SD")
repl_site   <- "OK"
disc_label  <- paste0("P < ", p_disc_cut)

setwd(readyData)
load("lipidomics_p4.Rda")
load("lipidomics_p5.Rda")
load("dat_AimOne_p4.Rda")
load("dat_AimOne_p5.Rda")

setwd(result)
P4P5all <- read.table('P4P5all.csv', sep=',', header=T)
UPF_lipids_all <- read.table('UPF_lipids_all.csv', sep=',', header=T)
lipids42 <- UPF_lipids_all$lipids_name

dat_AimOne_p4 <- dat_AimOne_p4 %>%
  select(-chd, -cvd, -TIME_chd, -TIME_cvd, -starts_with("incident"))

temp_1 <- as.data.frame(lipidomics_p4) %>%
  inner_join(dat_AimOne_p4, by = c("ID" = "IDNO"))
temp_2 <- as.data.frame(lipidomics_p5) %>%
  inner_join(dat_AimOne_p5, by = c("ID" = "IDNO"))

temp_1[,"Phase"] <- "Phase4"
temp_2[,"Phase"] <- "Phase5"
names(temp_1) <- sub('^S4', '', names(temp_1))
names(temp_2) <- sub('^S5', '', names(temp_2))
colnames(temp_1)[colnames(temp_1) == "PDM3"] <- "PA"
colnames(temp_2)[colnames(temp_2) == "ped_avg3day"] <- "PA"
temp_2$BOOKNUM <- NULL

temp_3 <- rbind(temp_1, temp_2)
temp_3$Phase  <- as.factor(temp_3$Phase)
temp_3$CENTER <- as.character(temp_3$CENTER)

first_ind <- which(colnames(temp_3) == "AC(10:0)")
last_ind  <- which(colnames(temp_3) == "Unknown_1230")
lipid_names <- colnames(temp_3)[first_ind:last_ind]

scale_vars <- c(grep("_P$", colnames(temp_3), value = TRUE), "AGE", "EDU", "BMI", "PA", "DT_KCAL")
scale_vars <- setdiff(scale_vars, lipid_names)

prep_set <- function(d){
  d %>%
    mutate(across(all_of(c(lipid_names, scale_vars)), ~ as.numeric(scale(.x, center = TRUE, scale = TRUE))),
           FID    = as.factor(as.character(FID)),
           ID     = as.factor(as.character(ID)),
           CENTER = as.factor(CENTER))
}

fit_all_lipids <- function(d, set_label){
  with_center <- length(unique(d$CENTER)) > 1
  rhs <- c("nova4_kcal_P", "AGE", "SEX", if(with_center) "CENTER",
           "EDU", "BMI", "SMOKE", "PA", "ETOH", "LIPIDLOW", "DT_KCAL", "Phase", "(1|FID)", "(1|ID)")
  fml <- as.formula(paste("lipid ~", paste(rhs, collapse = " + ")))
  d_cov <- d[, c("ID", "FID", "CENTER", "nova4_kcal_P", "AGE", "SEX", "EDU", "BMI", "SMOKE", "PA",
                 "ETOH", "LIPIDLOW", "DT_KCAL", "Phase")]
  out <- vector("list", length(lipid_names))
  for(i in seq_along(lipid_names)){
    d_cov$lipid <- d[[lipid_names[i]]]
    fit <- tryCatch(lmer(fml, data = d_cov), error = function(e) NULL)
    if(is.null(fit)){
      out[[i]] <- data.frame(lipids_name = lipid_names[i], beta = NA, se = NA, pvals = NA, N = NA)
      next
    }
    a <- coef(summary(fit))["nova4_kcal_P", ]
    out[[i]] <- data.frame(lipids_name = lipid_names[i], beta = a["Estimate"], se = a["Std. Error"],
                           pvals = a["Pr(>|t|)"], N = nobs(fit), row.names = NULL)
  }
  bind_rows(out) %>%
    mutate(lower_CI = beta - 1.96 * se,
           upper_CI = beta + 1.96 * se,
           qvals    = p.adjust(pvals, method = "BH"),
           known    = str_detect(lipids_name, "Unknown", negate = TRUE),
           set      = set_label)
}

setwd(result)
d_disc <- prep_set(temp_3 %>% filter(CENTER %in% disc_sites))
d_repl <- prep_set(temp_3 %>% filter(CENTER %in% repl_site))
disc <- fit_all_lipids(d_disc, paste0("discovery ", paste(disc_sites, collapse = "+")))
repl <- fit_all_lipids(d_repl, paste0("replication ", repl_site))
write.table(disc, "sensCenter_AZSD_OK_discovery_all.csv", sep=',', row.names = F, col.names=T)
write.table(repl, "sensCenter_AZSD_OK_replication_all.csv", sep=',', row.names = F, col.names=T)

full <- P4P5all %>% select(lipids_name, beta_full = beta, pvals_full = pvals, qvals_full = qvals)

comb <- disc %>%
  select(lipids_name, known, beta_disc = beta, se_disc = se, p_disc = pvals, N_disc = N) %>%
  left_join(repl %>% select(lipids_name, beta_repl = beta, se_repl = se, p_repl = pvals, N_repl = N),
            by = "lipids_name") %>%
  left_join(full, by = "lipids_name") %>%
  mutate(disc_center    = paste(disc_sites, collapse = "+"),
         repl_center    = repl_site,
         in42           = lipids_name %in% lipids42,
         disc_threshold = disc_label,
         disc_hit       = !is.na(p_disc) & p_disc < p_disc_cut,
         repl_primary   = coalesce(disc_hit & sign(beta_disc) == sign(beta_repl) & p_repl < p_repl_cut, FALSE)) %>%
  mutate(w_disc    = 1 / se_disc^2,
         w_repl    = 1 / se_repl^2,
         beta_meta = (w_disc * beta_disc + w_repl * beta_repl) / (w_disc + w_repl),
         se_meta   = sqrt(1 / (w_disc + w_repl)),
         lower_CI_meta = beta_meta - 1.96 * se_meta,
         upper_CI_meta = beta_meta + 1.96 * se_meta,
         p_meta    = 2 * pnorm(-abs(beta_meta / se_meta)),
         Q_het     = w_disc * (beta_disc - beta_meta)^2 + w_repl * (beta_repl - beta_meta)^2,
         p_het     = pchisq(Q_het, df = 1, lower.tail = FALSE),
         I2_het    = ifelse(Q_het > 1, (Q_het - 1) / Q_het, 0)) %>%
  select(-w_disc, -w_repl)
write.table(comb, "sensCenter_AZSD_OK_combined_all.csv", sep=',', row.names = F, col.names=T)

