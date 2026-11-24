# Aim 3: 5-year change in UPF-related lipids and incident MI and stroke after follow-up (Cox frailty).
# Input: 4.0 dataset. Output: result/ (used by 5.6, 8.1 and 8.1.1).

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
load("dat_lipidChange_P5landmark.Rda")
setwd(result)
key_lipids_risk <- read.table('UPF_lipids_risk.csv',sep=',',header=T)
key_lipids_prot <- read.table('UPF_lipids_prot.csv',sep=',',header=T)
lipid_sets <- list(riskLipids = key_lipids_risk$lipids_name,
                   protLipids = key_lipids_prot$lipids_name)

outcomes <- list(list(label = "Incident MI", short = "mi", atrisk = "atrisk_mi5", time_var = "TIME_mi5", event_var = "mi5"),
  list(label = "Incident definite stroke", short = "stroke", atrisk = "atrisk_stroke5", time_var = "TIME_stroke5", event_var = "stroke5"))

covs <- c("S5AGE","SEX","CENTER","S5EDU","S5BMI","S5SMOKE","S5ETOH","ped_avg3day",
          "LIPIDLOW_p5","S5ADADM2020","S5USHTN","S5DT_KCAL","BMI_diff","DT_KCAL_diff")
cont_vars <- c("S5AGE","S5EDU","S5BMI","ped_avg3day","S5DT_KCAL","BMI_diff","DT_KCAL_diff")

prep_sample <- function(out){
  d <- dat_lipidChange %>%
    filter(.data[[out$atrisk]]) %>%
    filter(if_all(all_of(c(covs, out$time_var)), ~ !is.na(.x)))
  for(v in cont_vars) d[[v]] <- as.numeric(scale(d[[v]]))
  sd_chg <- c()
  for(lp in unlist(lipid_sets)){
    sd_chg[lp] <- sd(d[[paste0(lp, "_chg")]], na.rm = TRUE)
    d[[paste0(lp, "_bsl")]] <- as.numeric(scale(d[[paste0(lp, "_bsl")]]))
    d[[paste0(lp, "_chg")]] <- as.numeric(scale(d[[paste0(lp, "_chg")]]))
  }
  d$FID <- as.factor(d$FID)
  list(data = d, sd_chg = sd_chg)
}

run_change_cox <- function(d, sd_chg, lipid_names, out, covs, model_label){
  n_lipids <- length(lipid_names)
  beta<-rep(NA,n_lipids); se<-rep(NA,n_lipids); HR<-rep(NA,n_lipids)
  lower_CI<-rep(NA,n_lipids); upper_CI<-rep(NA,n_lipids); pvals<-rep(NA,n_lipids)
  HR_bsl<-rep(NA,n_lipids); p_bsl<-rep(NA,n_lipids)
  form <- as.formula(paste0("Surv(", out$time_var, ",", out$event_var, ") ~ lipid_chg + lipid_bsl + ",
                            paste(covs, collapse = "+"), " + frailty(FID)"))
  for(i in 1:n_lipids){
    d$lipid_chg <- d[[paste0(lipid_names[i], "_chg")]]
    d$lipid_bsl <- d[[paste0(lipid_names[i], "_bsl")]]
    fit <- tryCatch(coxph(form, data = d),
                    error = function(e){  NULL })
    if(is.null(fit)) next
    a <- summary(fit)
    beta[i]     <- a$coefficients["lipid_chg", "coef"]
    se[i]       <- a$coefficients["lipid_chg", "se(coef)"]
    HR[i]       <- exp(beta[i])
    lower_CI[i] <- exp(beta[i] - 1.96 * se[i])
    upper_CI[i] <- exp(beta[i] + 1.96 * se[i])
    pvals[i]    <- a$coefficients["lipid_chg", "p"]
    HR_bsl[i]   <- exp(a$coefficients["lipid_bsl", "coef"])
    p_bsl[i]    <- a$coefficients["lipid_bsl", "p"]
  }
  qvals <- p.adjust(pvals, method = "BH")
  data.frame(lipids_name = lipid_names, outcome = out$label, model = model_label,
             N_total = nrow(d), N_events = sum(d[[out$event_var]]),
             beta, se, HR, lower_CI, upper_CI, pvals, qvals,
             HR_bsl, p_bsl,
             sd_change = as.numeric(sd_chg[lipid_names]))
}

setwd(result)
all_res <- list()
for(out in outcomes){
  s <- prep_sample(out)
  for(set_name in names(lipid_sets)){
    res <- run_change_cox(s$data, s$sd_chg, lipid_sets[[set_name]], out, covs, "adjusted")
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
    write.table(all_res[[keys[k]]], paste0('change_', out$short, '_', names(lipid_sets)[k], '.csv'),
                sep=',', row.names = F, col.names=T)
  }
}

res_all <- bind_rows(all_res)
write.table(res_all, 'change_mi_stroke_allLipids.csv', sep=',', row.names = F, col.names=T)
