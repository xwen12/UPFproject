# Center-split analysis: correlation of discovery and replication betas with the main betas, and power.
# Input: result/sensCenter_AZSD_OK_combined_all.csv (2.4). Output: result/sensCenter_OK_*.csv.

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
library(dplyr)

result <-file.path(proj_dir, "result")

setwd(result)
comb <- read.table('sensCenter_AZSD_OK_combined_all.csv', sep=',', header=T, check.names = FALSE)

kn <- comb %>% filter(known)

r_disc <- cor.test(kn$beta_disc, kn$beta_full, method = "pearson")
r_repl <- cor.test(kn$beta_repl, kn$beta_full, method = "pearson")

vs_main <- data.frame(
  half      = c("AZ+SD", "OK"),
  n_lipids  = nrow(kn),
  pearson_r = c(r_disc$estimate,    r_repl$estimate),
  lower_CI  = c(r_disc$conf.int[1], r_repl$conf.int[1]),
  upper_CI  = c(r_disc$conf.int[2], r_repl$conf.int[2]),
  p_value   = c(r_disc$p.value,     r_repl$p.value)) %>%
  mutate(across(c(pearson_r, lower_CI, upper_CI), ~ round(.x, 3)))

p42 <- kn %>%
  filter(in42) %>%
  mutate(z_disc     = abs(beta_full) / se_disc,
         z_repl     = abs(beta_full) / se_repl,
         power_disc = pnorm(z_disc - 1.96) + pnorm(-z_disc - 1.96),
         power_repl = pnorm(z_repl - 1.96),
         power_both = power_disc * power_repl,
         sig_disc   = p_disc < 0.05 & sign(beta_disc) == sign(beta_full),
         sig_repl   = p_repl < 0.05 & sign(beta_repl) == sign(beta_full),
         sig_both   = sig_disc & sig_repl)

pow_sum <- data.frame(
  measure    = c("P<0.05 in discovery (AZ+SD)", "P<0.05 in replication (OK)", "P<0.05 in both halves"),
  mean_power = c(mean(p42$power_disc), mean(p42$power_repl), mean(p42$power_both)),
  expected_n = c(sum(p42$power_disc),  sum(p42$power_repl),  sum(p42$power_both)),
  observed_n = c(sum(p42$sig_disc),    sum(p42$sig_repl),    sum(p42$sig_both))) %>%
  mutate(across(c(mean_power, expected_n), ~ round(.x, 3)))

setwd(result)
write.table(vs_main, 'sensCenter_OK_betaVsMain.csv',      sep=',', row.names = F, col.names=T)
write.table(pow_sum, 'sensCenter_OK_power42_summary.csv', sep=',', row.names = F, col.names=T)
