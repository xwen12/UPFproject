# Aim 1: baseline intake of each of 8 UPF subclasses and incident CVD, CHD, MI and stroke (Cox frailty).
# Input: analysis datasets (1.3) + 2.0.1 setup. Output: HR tables in result/.

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
raw_dir  <- "path/to/SHS_raw_data"
library(tidyverse)
library(survival)
rawData   <- raw_dir
readyData <- file.path(proj_dir, "data")
result    <- file.path(proj_dir, "result")
tabsFigs  <- file.path(proj_dir, "result", "tabsFigs")

setwd(readyData)
load("dat_AimOne_p4.Rda")
load("outcome_mi_stroke_p4.Rda")

temp_0 <- dat_AimOne_p4 %>%
  inner_join(outcome_mi_stroke_p4, by = "IDNO")

temp_1 <- temp_0 %>%
  mutate(
    nova4_1_10pct = nova4_1_kcal_P * 10,
    nova4_2_10pct = nova4_2_kcal_P * 10,
    nova4_3_10pct = nova4_3_kcal_P * 10,
    nova4_4_10pct = nova4_4_kcal_P * 10,
    nova4_5_10pct = nova4_5_kcal_P * 10,
    nova4_6_10pct = nova4_6_kcal_P * 10,
    nova4_7_10pct = nova4_7_kcal_P * 10,
    nova4_8_10pct = nova4_8_kcal_P * 10,
    S4AGE     = as.numeric(scale(S4AGE,     center = TRUE, scale = TRUE)),
    S4EDU     = as.numeric(scale(S4EDU,     center = TRUE, scale = TRUE)),
    S4BMI     = as.numeric(scale(S4BMI,     center = TRUE, scale = TRUE)),
    PDM3      = as.numeric(scale(PDM3,      center = TRUE, scale = TRUE)),
    S4DT_KCAL = as.numeric(scale(S4DT_KCAL, center = TRUE, scale = TRUE)),
    TIME_cvd  = as.numeric(TIME_cvd),
    TIME_chd  = as.numeric(TIME_chd),
    FID       = as.factor(FID)
  ) %>%
  arrange(FID)

base_covs <- c("S4AGE", "SEX", "CENTER", "S4EDU", "S4BMI",
               "S4SMOKE", "S4ETOH", "PDM3",
               "LIPIDLOW", "S4ADADM2020", "S4USHTN",
               "S4DT_KCAL",
               "frailty(FID)")

outcomes <- tibble(
  Outcome   = c("Incident CVD", "Incident CHD", "Incident MI", "Incident definite stroke"),
  time_var  = c("TIME_cvd",     "TIME_chd",     "TIME_mi",     "TIME_stroke"),
  event_var = c("cvd",          "chd",          "mi",          "definite_stroke"),
  Outcome_order = 1:4
)

exposures <- tibble(
  Exposure = paste0("nova4_", 1:8, "_10pct"),
  Exposure_label = c("Ultra-processed breads & breakfast foods",
                     "Sauces, spreads & condiments",
                     "Packaged sweet snacks & desserts",
                     "Packaged savory snacks",
                     "Sugar/artificially sweetened beverages",
                     "Meat/poultry/seafood-based ready-to-eat products",
                     "Ready-to-eat/heat mixed dishes",
                     "Dairy-based desserts"),
  Exposure_type  = "Subclass",
  Exposure_order = 1:8
)

fit_cox <- function(data, time_var, event_var, exposure_var){
  rhs  <- paste(c(exposure_var, base_covs), collapse = " + ")
  form <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ", rhs))
  fit <- tryCatch(coxph(form, data = data),
                  error = function(e){
                    NULL
                  })
  sm <- summary(fit)
  tibble(
    Exposure = exposure_var,
    N_used   = fit$n,
    N_events = fit$nevent,
    beta     = sm$coefficients[exposure_var, "coef"],
    se       = sm$coefficients[exposure_var, "se(coef)"],
    HR       = sm$conf.int[exposure_var, "exp(coef)"],
    Lower_CI = sm$conf.int[exposure_var, "lower .95"],
    Upper_CI = sm$conf.int[exposure_var, "upper .95"],
    P_value  = sm$coefficients[exposure_var, "p"]
  )
}

res_list <- list()
for (k in 1:nrow(outcomes)){
  for (j in 1:nrow(exposures)){
    res <- fit_cox(data         = temp_1,
                   time_var     = outcomes$time_var[k],
                   event_var    = outcomes$event_var[k],
                   exposure_var = exposures$Exposure[j])
    if (!is.null(res)) res_list[[paste(k, j)]] <- res %>% mutate(Outcome = outcomes$Outcome[k])
  }
}
results_raw <- bind_rows(res_list)

fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

results_subclass <- results_raw %>%
  group_by(Outcome) %>%
  mutate(n_tests = n(),
         q_value = p.adjust(P_value, method = "BH")) %>%
  ungroup() %>%
  left_join(outcomes %>% select(Outcome, Outcome_order), by = "Outcome") %>%
  left_join(exposures, by = "Exposure") %>%
  mutate(HR_95CI = paste0(sprintf("%.2f", HR), " (", sprintf("%.2f", Lower_CI),
                          ", ", sprintf("%.2f", Upper_CI), ")"),
         P_fmt   = fmt_p(P_value),
         q_fmt   = fmt_p(q_value),
         sig_q05 = q_value < 0.05) %>%
  arrange(Outcome_order, Exposure_order) %>%
  select(Outcome, Outcome_order, Exposure, Exposure_label, Exposure_type, Exposure_order,
         N_used, N_events, beta, se, HR, Lower_CI, Upper_CI,
         P_value, q_value, n_tests, HR_95CI, P_fmt, q_fmt, sig_q05)

setwd(result)
write.csv(results_subclass, "upf_subclass_cvd_outcomes.csv", row.names = FALSE)
