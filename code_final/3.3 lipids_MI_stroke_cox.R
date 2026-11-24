# Aim 3: baseline levels of UPF-related lipids and incident MI and definite stroke (Cox frailty models).
# Input: analysis datasets (1.3), 2.1 results. Output: result/ (used by 8.1 and 8.1.1).

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
raw_dir  <- "path/to/SHS_raw_data"
library(dplyr)
library(survival)
rawData <- raw_dir
readyData <- file.path(proj_dir, "data")
result <-file.path(proj_dir, "result")
tabsFigs <- file.path(proj_dir, "result", "tabsFigs")

setwd(readyData)
load("lipidomics_p4.Rda")
load("dat_AimOne_p4.Rda")
load("outcome_mi_stroke_p4.Rda")
setwd(result)
key_lipids_risk <- read.table('UPF_lipids_risk.csv',sep=',',header=T)
key_lipids_prot <- read.table('UPF_lipids_prot.csv',sep=',',header=T)

lipidomics_df <- as.data.frame(lipidomics_p4)
temp_1 <- lipidomics_df %>%
  inner_join(dat_AimOne_p4, by = c("ID" = "IDNO")) %>%
  inner_join(outcome_mi_stroke_p4, by = c("ID" = "IDNO"))

dt <- temp_1 %>%
  mutate(S4AGE     = as.numeric(scale(S4AGE)),
         S4EDU     = as.numeric(scale(S4EDU)),
         S4BMI     = as.numeric(scale(S4BMI)),
         PDM3      = as.numeric(scale(PDM3)),
         S4DT_KCAL = as.numeric(scale(S4DT_KCAL)),
         FID       = as.factor(FID))

all_key_lipids <- c(key_lipids_risk$lipids_name, key_lipids_prot$lipids_name)
for(lp in all_key_lipids){
  dt[[lp]] <- as.numeric(scale(dt[[lp]], center = TRUE, scale = TRUE))
}

dt <- dt %>% arrange(FID)

covs <- "S4AGE+SEX+CENTER+S4EDU+S4BMI+S4SMOKE+S4ETOH+PDM3+LIPIDLOW+S4ADADM2020+S4USHTN+S4DT_KCAL+frailty(FID)"

run_lipid_cox <- function(data, lipid_names, time_var, event_var, outcome_label){
  n_lipids <- length(lipid_names)
  beta<-rep(NA,n_lipids); se<-rep(NA,n_lipids); HR<-rep(NA,n_lipids)
  lower_CI<-rep(NA,n_lipids); upper_CI<-rep(NA,n_lipids); pvals<-rep(NA,n_lipids)
  form <- as.formula(paste0("Surv(", time_var, ",", event_var, ") ~ lipid+", covs))
  for(i in 1:n_lipids){
    data$lipid <- data[[lipid_names[i]]]
    fit <- tryCatch(coxph(form, data = data),
                    error = function(e){  NULL })
    if(is.null(fit)) next
    a <- summary(fit)
    beta[i]     <- a$coefficients["lipid", "coef"]
    se[i]       <- a$coefficients["lipid", "se(coef)"]
    HR[i]       <- exp(beta[i])
    lower_CI[i] <- exp(beta[i] - 1.96 * se[i])
    upper_CI[i] <- exp(beta[i] + 1.96 * se[i])
    pvals[i]    <- a$coefficients["lipid", "p"]
  }
  qvals <- p.adjust(pvals, method = "BH")
  data.frame(lipids_name = lipid_names,
             outcome     = outcome_label,
             N_total     = nrow(data),
             N_events    = sum(data[[event_var]] == 1),
             beta, se, HR, lower_CI, upper_CI, pvals, qvals)
}

outcomes <- list(
  list(label = "Incident MI",              short = "mi",     time_var = "TIME_mi",     event_var = "mi"),
  list(label = "Incident definite stroke", short = "stroke", time_var = "TIME_stroke", event_var = "definite_stroke"))
lipid_sets <- list(riskLipids = key_lipids_risk$lipids_name,
                   protLipids = key_lipids_prot$lipids_name)

setwd(result)
all_res <- list()
for(out in outcomes){
  for(set_name in names(lipid_sets)){
    res <- run_lipid_cox(dt, lipid_sets[[set_name]], out$time_var, out$event_var, out$label)
    res$lipid_set <- set_name
    all_res[[paste(out$short, set_name)]] <- res
  }
  keys  <- paste(out$short, names(lipid_sets))
  q_all <- p.adjust(unlist(lapply(all_res[keys], `[[`, "pvals")), method = "BH")
  n_set <- sapply(all_res[keys], nrow)
  q_split <- split(q_all, rep(seq_along(keys), n_set))
  for(k in seq_along(keys)){
    all_res[[keys[k]]]$qvals_withinSet <- all_res[[keys[k]]]$qvals
    all_res[[keys[k]]]$qvals <- q_split[[k]]
    write.table(all_res[[keys[k]]], paste0('longitu_', out$short, '_', names(lipid_sets)[k], '.csv'),
                sep=',', row.names = F, col.names=T)
  }
}

res_all <- bind_rows(all_res)
write.table(res_all,'longitu_mi_stroke_allLipids.csv',sep=',',row.names = F, col.names=T)
