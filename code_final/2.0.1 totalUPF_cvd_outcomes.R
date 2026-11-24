# Aim 1: baseline total UPF intake and incident CVD, CHD, MI and definite stroke (Cox frailty models).
# Input: analysis datasets (1.3). Output: HR tables in result/.

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
raw_dir  <- "path/to/SHS_raw_data"
library(tidyverse)
library(survival)
library(lubridate)
rawData   <- raw_dir
readyData <- file.path(proj_dir, "data")
result    <- file.path(proj_dir, "result")
tabsFigs  <- file.path(proj_dir, "result", "tabsFigs")

setwd(readyData)
load("dat_AimOne_p4.Rda")

setwd(rawData)
chd_raw <- read.csv("cvd data_2020.csv", header = TRUE, stringsAsFactors = FALSE)

closing_date    <- as.Date("2020-12-31")

parse_shs_date <- function(x){
  x <- as.character(x)
  x[x == ""] <- NA
  as.Date(parse_date_time(x, orders = c("mdY", "Ymd", "dbY"), quiet = TRUE))
}

cvd_subtypes <- chd_raw %>%
  select(IDNO, S4EXAMDT, ANYMIS4, ANYMIDTS4, ANYDSTKS4, ANYDSTKDTS4, DOD) %>%
  mutate(ANYMIS4     = as.integer(ANYMIS4),
         ANYDSTKS4   = as.integer(ANYDSTKS4),
         S4EXAMDT    = parse_shs_date(S4EXAMDT),
         ANYMIDTS4   = parse_shs_date(ANYMIDTS4),
         ANYDSTKDTS4 = parse_shs_date(ANYDSTKDTS4),
         DOD         = parse_shs_date(DOD))

temp_0 <- dat_AimOne_p4 %>%
  left_join(cvd_subtypes, by = "IDNO")

temp_0 <- temp_0 %>%
  mutate(
    end_fu = pmin(DOD, closing_date, na.rm = TRUE),
    mi      = ifelse(!is.na(ANYMIS4) & ANYMIS4 == 1, 1L, 0L),
    mi_date = case_when(mi == 1 & !is.na(ANYMIDTS4) ~ ANYMIDTS4,
                        mi == 1 &  is.na(ANYMIDTS4) ~ DOD,
                        TRUE                        ~ end_fu),
    TIME_mi = as.numeric(mi_date - S4EXAMDT),
    definite_stroke = ifelse(!is.na(ANYDSTKS4) & ANYDSTKS4 == 1, 1L, 0L),
    stroke_date     = case_when(definite_stroke == 1 & !is.na(ANYDSTKDTS4) ~ ANYDSTKDTS4,
                                definite_stroke == 1 &  is.na(ANYDSTKDTS4) ~ DOD,
                                TRUE                                       ~ end_fu),
    TIME_stroke = as.numeric(stroke_date - S4EXAMDT)
  )

temp_0 <- temp_0 %>%
  filter(!is.na(TIME_mi), !is.na(TIME_stroke), TIME_mi > 0, TIME_stroke > 0)

outcome_mi_stroke_p4 <- temp_0 %>%
  select(IDNO, mi, TIME_mi, definite_stroke, TIME_stroke)
setwd(readyData)
save(outcome_mi_stroke_p4, file = "outcome_mi_stroke_p4.Rda")

temp_1 <- temp_0 %>%
  mutate(
    nova4_10pct = nova4_kcal_P * 10,
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
  res <- fit_cox(data         = temp_1,
                 time_var     = outcomes$time_var[k],
                 event_var    = outcomes$event_var[k],
                 exposure_var = "nova4_10pct")
  if (!is.null(res)) res_list[[k]] <- res %>% mutate(Outcome = outcomes$Outcome[k])
}
results_raw <- bind_rows(res_list)

fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

results_total <- results_raw %>%
  mutate(n_tests = n(),
         q_value = p.adjust(P_value, method = "BH")) %>%
  left_join(outcomes %>% select(Outcome, Outcome_order), by = "Outcome") %>%
  mutate(Exposure_label = "Total UPF",
         Exposure_type  = "Total",
         HR_95CI = paste0(sprintf("%.2f", HR), " (", sprintf("%.2f", Lower_CI),
                          ", ", sprintf("%.2f", Upper_CI), ")"),
         P_fmt   = fmt_p(P_value),
         q_fmt   = fmt_p(q_value),
         sig_q05 = q_value < 0.05) %>%
  arrange(Outcome_order) %>%
  select(Outcome, Outcome_order, Exposure, Exposure_label, Exposure_type,
         N_used, N_events, beta, se, HR, Lower_CI, Upper_CI,
         P_value, q_value, n_tests, HR_95CI, P_fmt, q_fmt, sig_q05)

setwd(result)
write.csv(results_total, "upf_total_cvd_outcomes.csv", row.names = FALSE)
