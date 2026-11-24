# Supplemental Table S2: total UPF and incident CVD outcomes, main model vs AHEI-adjusted model.
# Input: 2.0.1 and 2.2 results. Output: Table S2 (result/tabsFigs).

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
library(tidyverse)
library(openxlsx)
result    <- file.path(proj_dir, "result")
tabsFigs  <- file.path(proj_dir, "result", "tabsFigs")

if (!dir.exists(tabsFigs)) dir.create(tabsFigs, recursive = TRUE)

setwd(result)
main_raw <- read.csv("upf_total_cvd_outcomes.csv", stringsAsFactors = FALSE)
ahei_raw <- read.csv("sensitivity_upf_ahei_long.csv", stringsAsFactors = FALSE)

upf_MB <- ahei_raw %>% filter(Term == "nova4_10pct", Model == "Model B: + AHEI")

outcome_order <- c("Incident CVD", "Incident CHD", "Incident MI", "Incident definite stroke")
outcome_names <- c("Incident CVD"             = "CVD",
                   "Incident CHD"             = "CHD",
                   "Incident MI"              = "MI",
                   "Incident definite stroke" = "Definite stroke")

fmt_p  <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
fmt_hr <- function(hr, lo, hi) sprintf("%.2f (%.2f, %.2f)", hr, lo, hi)
fmt_n  <- function(x) format(x, big.mark = ",", trim = TRUE)

note_hr    <- paste0("Hazard ratios (HR) and 95% confidence intervals (CI) are per 10% of total energy intake from UPF, ",
                     "estimated using Cox proportional hazards models with a shared frailty term for family.")
note_main  <- paste0("Main model: adjusted for age, sex, study center, education, body mass index, smoking, ",
                     "alcohol intake, physical activity, lipid-lowering medication use, glycemic status, ",
                     "hypertension status and total energy intake.")
note_ahei  <- paste0("AHEI-adjusted model: main model further adjusted for the Alternative Healthy Eating ",
                     "Index-2010 (AHEI, per SD), among participants with available AHEI data.")
note_fdr   <- "FDR: Benjamini-Hochberg false discovery rate across the 4 outcomes, calculated separately for each model."
note_abbr  <- paste0("Abbreviations: AHEI, Alternative Healthy Eating Index; CHD, coronary heart disease; ",
                     "CVD, cardiovascular disease; MI, myocardial infarction; UPF, ultra-processed food.")
notes_S2   <- c(note_hr, note_main, note_ahei, note_fdr, note_abbr)

main_num <- main_raw %>%
  transmute(Outcome, Model = "Main",
            N_used, N_events, beta, HR, Lower_CI, Upper_CI, P_value, q_value)

ahei_num <- upf_MB %>%
  mutate(q_value = p.adjust(P_value, method = "BH")) %>%
  transmute(Outcome, Model = "AHEI",
            N_used, N_events, beta, HR, Lower_CI, Upper_CI, P_value, q_value)

s2_num <- bind_rows(main_num, ahei_num) %>%
  mutate(Outcome = factor(Outcome, levels = outcome_order),
         Model   = factor(Model,   levels = c("Main", "AHEI"))) %>%
  arrange(Model, Outcome)

tabS2 <- s2_num %>%
  transmute(Outcome,
            Model,
            cases_n = paste0(fmt_n(N_events), " / ", fmt_n(N_used)),
            hr_ci   = fmt_hr(HR, Lower_CI, Upper_CI),
            p       = fmt_p(P_value),
            fdr     = fmt_p(q_value)) %>%
  pivot_wider(names_from = Model, values_from = c(cases_n, hr_ci, p, fdr)) %>%
  arrange(Outcome) %>%
  transmute(Outcome       = unname(outcome_names[as.character(Outcome)]),
            cases_n_Main, hr_ci_Main, p_Main, fdr_Main,
            cases_n_AHEI, hr_ci_AHEI, p_AHEI, fdr_AHEI)

header_top <- c("", "Main model", "", "", "", "AHEI-adjusted model", "", "", "")
header_sub <- c("Outcome", "Cases / N", "HR (95% CI)", "P", "FDR",
                           "Cases / N", "HR (95% CI)", "P", "FDR")

wb <- createWorkbook()
addWorksheet(wb, "Table S2")
n_col <- ncol(tabS2)

writeData(wb, "Table S2", t(header_top), startRow = 1, startCol = 1, colNames = FALSE)
writeData(wb, "Table S2", t(header_sub), startRow = 2, startCol = 1, colNames = FALSE)
writeData(wb, "Table S2", tabS2,         startRow = 3, startCol = 1, colNames = FALSE)

mergeCells(wb, "Table S2", cols = 2:5, rows = 1)
mergeCells(wb, "Table S2", cols = 6:9, rows = 1)

body_end   <- 2 + nrow(tabS2)
note_start <- body_end + 2
writeData(wb, "Table S2", notes_S2, startRow = note_start, startCol = 1)

top_style  <- createStyle(textDecoration = "bold", halign = "center", valign = "center",
                          border = "bottom")
sub_style  <- createStyle(textDecoration = "bold", halign = "center", valign = "center",
                          border = "bottom", wrapText = TRUE)
body_style <- createStyle(halign = "center", valign = "center")
out_style  <- createStyle(halign = "left",   valign = "center")
note_style <- createStyle(wrapText = TRUE, valign = "top", fontSize = 10)
addStyle(wb, "Table S2", top_style,  rows = 1, cols = 2:n_col, gridExpand = TRUE)
addStyle(wb, "Table S2", sub_style,  rows = 2, cols = 1:n_col, gridExpand = TRUE)
addStyle(wb, "Table S2", body_style, rows = 3:body_end, cols = 2:n_col, gridExpand = TRUE)
addStyle(wb, "Table S2", out_style,  rows = 3:body_end, cols = 1, gridExpand = TRUE)
addStyle(wb, "Table S2", createStyle(border = "bottom"), rows = body_end, cols = 1:n_col,
         gridExpand = TRUE, stack = TRUE)

for (r in note_start:(note_start + length(notes_S2) - 1)) {
  mergeCells(wb, "Table S2", cols = 1:n_col, rows = r)
  addStyle(wb, "Table S2", note_style, rows = r, cols = 1:n_col, gridExpand = TRUE)
}
setRowHeights(wb, "Table S2", rows = note_start:(note_start + length(notes_S2) - 1), heights = 30)

setColWidths(wb, "Table S2", cols = 1:n_col, widths = c(16, 12, 18, 8, 8, 12, 18, 8, 8))

setwd(tabsFigs)
saveWorkbook(wb, "tableS2_totalUPF_AHEI.xlsx", overwrite = TRUE)
setwd(result)
write.csv(s2_num, "tableS2_totalUPF_AHEI_numbers.csv", row.names = FALSE)
