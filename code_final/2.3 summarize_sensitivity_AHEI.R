# Compare hazard ratios between the main and AHEI-adjusted models (2.2).
# Input: result/sensitivity_upf_ahei_*.csv. Output: result/sensitivity_upf_ahei_compare.csv.

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
library(tidyverse)

result <- file.path(proj_dir, "result")

setwd(result)
res_long  <- read.csv("sensitivity_upf_ahei_long.csv",  stringsAsFactors = FALSE, check.names = FALSE)
res_table <- read.csv("sensitivity_upf_ahei_table.csv", stringsAsFactors = FALSE, check.names = FALSE)

outcome_order <- c("Incident CVD", "Incident CHD", "Incident MI", "Incident definite stroke")

upf <- res_long %>%
  filter(Term == "nova4_10pct") %>%
  mutate(Model = recode(Model,
                        "Model 0: main, full sample"     = "M0",
                        "Model A: main, AHEI sub-sample" = "MA",
                        "Model B: + AHEI"                = "MB")) %>%
  select(Outcome, Model, N_used, N_events, beta, HR, Lower_CI, Upper_CI, P_value, HR_95CI)

upf_wide <- upf %>%
  pivot_wider(names_from = Model,
              values_from = c(N_used, N_events, beta, HR, Lower_CI, Upper_CI, P_value, HR_95CI)) %>%
  mutate(Outcome = factor(Outcome, levels = outcome_order)) %>%
  arrange(Outcome)

compare <- upf_wide %>%
  mutate(
    sig_M0         = P_value_M0 < 0.05,
    sig_MB         = P_value_MB < 0.05,
    same_direction = sign(beta_M0) == sign(beta_MB),
    kept           = sig_M0 & sig_MB & same_direction,
    new_sig        = !sig_M0 & sig_MB,
    pct_0toB       = (beta_MB - beta_M0) / beta_M0 * 100,
    pct_0toA       = (beta_MA - beta_M0) / beta_M0 * 100,
    pct_AtoB       = (beta_MB - beta_MA) / beta_MA * 100,
    beta_near_0    = abs(beta_M0) < 0.02,
    within_10pct   = abs(pct_AtoB) <= 10
  )

write.csv(compare, "sensitivity_upf_ahei_compare.csv", row.names = FALSE)
