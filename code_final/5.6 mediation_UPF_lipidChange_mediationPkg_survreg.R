# Aim 3 mediation: UPF -> 5-year lipid change -> incident outcomes (mediation package, Weibull survreg).
# Candidate mediators: lipid changes with P < 0.05 in 4.1-4.3. Main and sensitivity (no BMI, glycemic status, hypertension) models.
# Output: result/mediation_change_mediationPkg_survreg_*.csv.

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
raw_dir  <- "path/to/SHS_raw_data"
library(mediation)
library(sandwich)
library(survival)
library(dplyr)
rawData <- raw_dir
readyData <- file.path(proj_dir, "data")
result <-file.path(proj_dir, "result")
tabsFigs <- file.path(proj_dir, "result", "tabsFigs")

sig_cut      <- 0.05
n_sims       <- 1000
surv_dist    <- "weibull"
time_unit    <- 365.25
run_covset   <- c("main", "sensNoBMIDMHTN")

covset_list <- list(
  main           = list(tag = "_main",
                        vars = c("S4AGE","SEX","CENTER","S4EDU","S4BMI","S4SMOKE","S4ETOH","PDM3",
                                 "LIPIDLOW","S4ADADM2020","S4USHTN","S4DT_KCAL")),
  sensNoBMIDMHTN = list(tag = "_sensNoBMIDMHTN",
                        vars = c("S4AGE","SEX","CENTER","S4EDU","S4SMOKE","S4ETOH","PDM3",
                                 "LIPIDLOW","S4DT_KCAL")))
sample_vars <- covset_list$main$vars
cont_vars   <- c("S4AGE","S4EDU","S4BMI","PDM3","S4DT_KCAL")

setwd(readyData)
load("dat_AimOne_p4.Rda")
load("dat_lipidChange_P5landmark.Rda")

out_chg <- list(
  list(label = "Incident CVD", short = "cvd", atrisk = "atrisk_cvd5", time_var = "TIME_cvd5", event_var = "cvd5",
       res_files = c("change_cvd_riskLipids.csv", "change_cvd_protLipids.csv")),
  list(label = "Incident CHD", short = "chd", atrisk = "atrisk_chd5", time_var = "TIME_chd5", event_var = "chd5",
       res_files = c("change_chd_riskLipids.csv", "change_chd_protLipids.csv")),
  list(label = "Incident MI",  short = "mi",  atrisk = "atrisk_mi5",  time_var = "TIME_mi5",  event_var = "mi5",
       res_files = c("change_mi_riskLipids.csv",  "change_mi_protLipids.csv")))

setwd(result)
for(k in seq_along(out_chg)){
  r <- bind_rows(lapply(out_chg[[k]]$res_files, read.csv))
  out_chg[[k]]$lipids <- unique(r$lipids_name[!is.na(r$pvals) & r$pvals < sig_cut])
}

p4_vars <- dat_AimOne_p4 %>%
  dplyr::select(IDNO, nova4_kcal_P, S4AGE, S4EDU, S4BMI, S4SMOKE, S4ETOH, PDM3,
                LIPIDLOW, S4ADADM2020, S4USHTN, S4DT_KCAL)
dat_chg <- dat_lipidChange %>%
  inner_join(p4_vars, by = "IDNO") %>%
  mutate(nova4_10pct = nova4_kcal_P * 10)

prep_sample <- function(out){
  med_cols <- c(paste0(out$lipids, "_chg"), paste0(out$lipids, "_bsl"))
  d <- dat_chg %>%
    filter(.data[[out$atrisk]]) %>%
    filter(if_all(all_of(c("nova4_10pct", sample_vars, out$time_var, med_cols)), ~ !is.na(.x)))
  for(v in cont_vars) d[[v]] <- as.numeric(scale(d[[v]]))
  for(v in med_cols)  d[[v]] <- as.numeric(scale(d[[v]]))
  d$event   <- as.integer(d[[out$event_var]])
  d$time_yr <- as.numeric(d[[out$time_var]]) / time_unit
  n_nonpos  <- sum(d$time_yr <= 0)
  if(n_nonpos > 0){
    d$time_yr[d$time_yr <= 0] <- 1 / time_unit
  }
  d$FID <- as.factor(d$FID)
  d
}

run_one <- function(d, lipid_name, out, cov_vars){
  covs <- paste(cov_vars, collapse = "+")

  d$lipid     <- d[[paste0(lipid_name, "_chg")]]
  d$lipid_bsl <- d[[paste0(lipid_name, "_bsl")]]

  d <- d[complete.cases(d[, c("lipid", "lipid_bsl", "event", "time_yr", "nova4_10pct", sample_vars)]), ]
  rownames(d) <- NULL

  f_m <- as.formula(paste0("lipid ~ nova4_10pct+lipid_bsl+", covs))
  f_y <- as.formula(paste0("Surv(time_yr, event) ~ nova4_10pct+lipid+lipid_bsl+", covs, "+cluster(FID)"))
  fit_m <- eval(bquote(lm(.(f_m), data = d)))
  fit_y <- eval(bquote(survreg(.(f_y), data = d, dist = .(surv_dist))))

  ctrl <- median(d$nova4_10pct)
  trt  <- ctrl + 1

  set.seed(1234)
  med <- tryCatch(
    withCallingHandlers(
      mediate(fit_m, fit_y, treat = "nova4_10pct", mediator = "lipid",
              control.value = ctrl, treat.value = trt, sims = n_sims, cluster = d$FID),
      warning = function(w){
        if(grepl("cluster' ignored for survival models", conditionMessage(w))) invokeRestart("muffleWarning")
      }),
    error = function(e){  NULL })

  v_m   <- sandwich::vcovCL(fit_m, cluster = d$FID)
  a_est <- unname(coef(fit_m)["nova4_10pct"]); a_se <- sqrt(v_m["nova4_10pct", "nova4_10pct"])
  b_est <- unname(coef(fit_y)["lipid"]);       b_se <- sqrt(vcov(fit_y)["lipid", "lipid"])
  c_est <- unname(coef(fit_y)["nova4_10pct"]); c_se <- sqrt(vcov(fit_y)["nova4_10pct", "nova4_10pct"])
  data.frame(
    analysis = "change", outcome = out$label, lipids_name = lipid_name,
    N_total = nrow(d), N_events = sum(d$event), fu_years_median = median(d$time_yr),
    upf_contrast = paste0(round(ctrl, 2), " -> ", round(trt, 2)),
    a_beta = a_est, a_lo = a_est - 1.96 * a_se, a_hi = a_est + 1.96 * a_se, a_p = 2 * pnorm(-abs(a_est / a_se)),
    b_TR   = exp(b_est), b_TR_lo = exp(b_est - 1.96 * b_se), b_TR_hi = exp(b_est + 1.96 * b_se),
    b_p    = 2 * pnorm(-abs(b_est / b_se)),
    c_TR   = exp(c_est), c_TR_lo = exp(c_est - 1.96 * c_se), c_TR_hi = exp(c_est + 1.96 * c_se),
    ACME = med$d.avg,    ACME_lo = med$d.avg.ci[1], ACME_hi = med$d.avg.ci[2], ACME_p = med$d.avg.p,
    ADE  = med$z.avg,    ADE_lo  = med$z.avg.ci[1], ADE_hi  = med$z.avg.ci[2], ADE_p  = med$z.avg.p,
    TE   = med$tau.coef, TE_lo   = med$tau.ci[1],   TE_hi   = med$tau.ci[2],   TE_p   = med$tau.p,
    PropMed = med$n.avg, PropMed_lo = med$n.avg.ci[1], PropMed_hi = med$n.avg.ci[2],
    same_sign = sign(med$d.avg) == sign(med$z.avg),
    dist = surv_dist, weibull_scale = fit_y$scale, covset = cs_name,
    covariates = covs, row.names = NULL)
}

for(cs_name in run_covset){
  cov_vars <- covset_list[[cs_name]]$vars

  res_list <- list(); fit_check <- list()
  for(out in out_chg){
    lipid_names <- out$lipids
    if(length(lipid_names) == 0){  next }

    d_out <- prep_sample(out)

    f_base <- as.formula(paste0("Surv(time_yr, event) ~ nova4_10pct+", paste(cov_vars, collapse = "+")))
    aic <- sapply(c("weibull", "lognormal"), function(ds)
      AIC(eval(bquote(survreg(.(f_base), data = d_out, dist = .(ds))))))
    fit_check[[out$short]] <- data.frame(outcome = out$label, AIC_weibull = aic[["weibull"]],
                                         AIC_lognormal = aic[["lognormal"]])

    res <- NULL
    for(i in seq_along(lipid_names)){
      res <- rbind(res, run_one(d_out, lipid_names[i], out, cov_vars))
    }
    if(is.null(res)) next

    res$ACME_q <- p.adjust(res$ACME_p, method = "BH")
    res_list[[out$short]] <- res
  }
  res_all <- bind_rows(res_list)
  fit_all <- bind_rows(fit_check)

  summ <- res_all %>%
    group_by(analysis, outcome) %>%
    summarise(n_tested       = n(),
              n_ACME_p05     = sum(ACME_p < 0.05),
              n_ACME_q05     = sum(ACME_q < 0.05),
              n_p05_samesign = sum(ACME_p < 0.05 & same_sign),
              TE_p           = round(first(TE_p), 3),
              .groups = "drop") %>%
    left_join(fit_all, by = "outcome")

  setwd(result)
  tag <- paste0("_mediationPkg_survreg", covset_list[[cs_name]]$tag)
  write.table(res_all, paste0("mediation_change", tag, ".csv"), sep = ",", row.names = F, col.names = T)
  write.table(summ,    paste0("mediation_change_summary", tag, ".csv"), sep = ",", row.names = F, col.names = T)
}
