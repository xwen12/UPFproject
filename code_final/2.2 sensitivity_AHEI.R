# Sensitivity analysis: total UPF and incident CVD outcomes, further adjusted for diet quality (AHEI).
# Input: analysis datasets (1.3) + AHEI. Output: result/sensitivity_upf_ahei_*.csv.

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

setwd(rawData)
ahei <- read.csv("p4nutrients/p4ahei.csv", header = TRUE, stringsAsFactors = FALSE)
ahei <- ahei %>% select(IDNO, AHEI)

temp_0 <- dat_AimOne_p4 %>%
  inner_join(outcome_mi_stroke_p4, by = "IDNO")

temp_0 <- temp_0 %>%
  left_join(ahei, by = "IDNO") %>%
  mutate(has_ahei = !is.na(AHEI))

temp_1 <- temp_0 %>%
  mutate(
    nova4_10pct = nova4_kcal_P * 10,
    S4AGE     = as.numeric(scale(S4AGE,     center = TRUE, scale = TRUE)),
    S4EDU     = as.numeric(scale(S4EDU,     center = TRUE, scale = TRUE)),
    S4BMI     = as.numeric(scale(S4BMI,     center = TRUE, scale = TRUE)),
    PDM3      = as.numeric(scale(PDM3,      center = TRUE, scale = TRUE)),
    S4DT_KCAL = as.numeric(scale(S4DT_KCAL, center = TRUE, scale = TRUE)),
    AHEI_z    = as.numeric(scale(AHEI,      center = TRUE, scale = TRUE)),
    FID       = as.factor(FID)
  ) %>%
  arrange(FID)

temp_ahei <- temp_1 %>% filter(has_ahei)

base_covs <- c("S4AGE", "SEX", "CENTER", "S4EDU", "S4BMI",
               "S4SMOKE", "S4ETOH", "PDM3",
               "LIPIDLOW", "S4ADADM2020", "S4USHTN",
               "S4DT_KCAL",
               "frailty(FID)")

outcomes <- list(
  list(label = "Incident CVD",             time_var = "TIME_cvd",    event_var = "cvd"),
  list(label = "Incident CHD",             time_var = "TIME_chd",    event_var = "chd"),
  list(label = "Incident MI",              time_var = "TIME_mi",     event_var = "mi"),
  list(label = "Incident definite stroke", time_var = "TIME_stroke", event_var = "definite_stroke")
)

fit_cox <- function(data, time_var, event_var, extra_covs = NULL,
                    outcome_label, model_label) {

  rhs  <- paste(c("nova4_10pct", extra_covs, base_covs), collapse = " + ")
  form <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ", rhs))

  fit <- tryCatch(coxph(form, data = data),
                  error = function(e) {
                    NULL
                  })

  sm <- summary(fit)
  terms <- intersect(c("nova4_10pct", "AHEI_z"), rownames(sm$conf.int))

  data.frame(
    Outcome  = outcome_label,
    Model    = model_label,
    Term     = terms,
    N_used   = fit$n,
    N_events = fit$nevent,
    beta     = sm$coefficients[terms, "coef"],
    HR       = sm$conf.int[terms, "exp(coef)"],
    Lower_CI = sm$conf.int[terms, "lower .95"],
    Upper_CI = sm$conf.int[terms, "upper .95"],
    P_value  = sm$coefficients[terms, "p"],
    row.names = NULL
  )
}

res_list <- list()
for (out in outcomes) {
  res_list[[paste0(out$label, "_0")]] <- fit_cox(temp_1, out$time_var, out$event_var,
                                                 extra_covs = NULL, out$label,
                                                 "Model 0: main")
  res_list[[paste0(out$label, "_B")]] <- fit_cox(temp_ahei, out$time_var, out$event_var,
                                                 extra_covs = "AHEI_z", out$label,
                                                 "Model B: + AHEI")
}

results_long <- bind_rows(res_list) %>%
  mutate(HR_95CI = paste0(sprintf("%.2f", HR), " (", sprintf("%.2f", Lower_CI),
                          "\u2013", sprintf("%.2f", Upper_CI), ")"))

pct_change <- results_long %>%
  filter(Term == "nova4_10pct") %>%
  select(Outcome, Model, beta) %>%
  pivot_wider(names_from = Model, values_from = beta) %>%
  rename(beta_0 = `Model 0: main`, beta_B = `Model B: + AHEI`) %>%
  mutate(pct_change = round((beta_B - beta_0) / beta_0 * 100, 1))

q_B <- results_long %>%
  filter(Term == "nova4_10pct", Model == "Model B: + AHEI") %>%
  mutate(q_value_B = p.adjust(P_value, method = "BH")) %>%
  select(Outcome, q_value_B)

fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

upf_wide <- results_long %>%
  filter(Term == "nova4_10pct") %>%
  mutate(Model = recode(Model, "Model 0: main" = "M0", "Model B: + AHEI" = "MB"),
         P_fmt = fmt_p(P_value)) %>%
  select(Outcome, Model, N_used, N_events, HR_95CI, P_fmt) %>%
  pivot_wider(names_from = Model, values_from = c(N_used, N_events, HR_95CI, P_fmt))

ahei_hr <- results_long %>%
  filter(Term == "AHEI_z") %>%
  transmute(Outcome, `AHEI HR per SD (95% CI)` = HR_95CI, `AHEI P` = fmt_p(P_value))

table_ahei <- upf_wide %>%
  left_join(pct_change %>% select(Outcome, pct_change), by = "Outcome") %>%
  left_join(q_B,     by = "Outcome") %>%
  left_join(ahei_hr, by = "Outcome") %>%
  transmute(
    Outcome,
    `Model 0 N/cases`         = paste0(N_used_M0, "/", N_events_M0),
    `Model 0 HR (95% CI)`     = HR_95CI_M0,
    `Model 0 P`               = P_fmt_M0,
    `Model B N/cases`         = paste0(N_used_MB, "/", N_events_MB),
    `Model B HR (95% CI)`     = HR_95CI_MB,
    `Model B P`               = P_fmt_MB,
    `Model B FDR q`           = sprintf("%.3f", q_value_B),
    `% change log-HR (0 to B)` = pct_change,
    `AHEI HR per SD (95% CI)`,
    `AHEI P`
  ) %>%
  mutate(Outcome = factor(Outcome, levels = sapply(outcomes, `[[`, "label"))) %>%
  arrange(Outcome)

setwd(result)
write.csv(table_ahei,   "sensitivity_upf_ahei_table.csv", row.names = FALSE)
write.csv(results_long, "sensitivity_upf_ahei_long.csv",  row.names = FALSE)
