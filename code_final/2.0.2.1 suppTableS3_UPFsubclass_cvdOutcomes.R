# Supplemental Table S3: UPF subclasses and incident CVD, CHD, MI and stroke.
# Input: 2.0.2 results. Output: Table S3 (result/tabsFigs).

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
library(tidyverse)
library(openxlsx)
result    <- file.path(proj_dir, "result")
tabsFigs  <- file.path(proj_dir, "result", "tabsFigs")

if (!dir.exists(tabsFigs)) dir.create(tabsFigs, recursive = TRUE)

setwd(result)
sub_raw <- read.csv("upf_subclass_cvd_outcomes.csv", stringsAsFactors = FALSE)

outcome_order <- c("Incident CVD", "Incident CHD", "Incident MI", "Incident definite stroke")
outcome_names <- c("Incident CVD"             = "CVD",
                   "Incident CHD"             = "CHD",
                   "Incident MI"              = "MI",
                   "Incident definite stroke" = "Definite stroke")
outcome_keys  <- c("Incident CVD"             = "cvd",
                   "Incident CHD"             = "chd",
                   "Incident MI"              = "mi",
                   "Incident definite stroke" = "stroke")

cases_check <- sub_raw %>%
  group_by(Outcome) %>%
  summarise(n_diff_events = n_distinct(N_events),
            n_diff_N      = n_distinct(N_used),
            min_events = min(N_events), max_events = max(N_events),
            min_N      = min(N_used),   max_N      = max(N_used),
            .groups = "drop") %>%
  mutate(same_cases = n_diff_events == 1 & n_diff_N == 1)
all_same_cases <- all(cases_check$same_cases)

fmt_p  <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
fmt_hr <- function(hr, lo, hi) sprintf("%.2f (%.2f, %.2f)", hr, lo, hi)
fmt_n  <- function(x) format(x, big.mark = ",", trim = TRUE)

cases_lab <- cases_check %>%
  mutate(Outcome = factor(Outcome, levels = outcome_order),
         label   = ifelse(same_cases,
                          paste0("Cases / N: ", fmt_n(min_events), " / ", fmt_n(min_N)),
                          paste0("Cases / N: ", fmt_n(min_events), "-", fmt_n(max_events),
                                 " / ", fmt_n(min_N), "-", fmt_n(max_N)))) %>%
  arrange(Outcome)

note_hr    <- paste0("Hazard ratios (HR) and 95% confidence intervals (CI) are per 10% of total energy intake ",
                     "from each UPF subclass, estimated using Cox proportional hazards models with a shared ",
                     "frailty term for family. Each subclass was analyzed in a separate model (subclasses were ",
                     "not mutually adjusted).")
note_main  <- paste0("All models were adjusted for age, sex, study center, education, body mass index, smoking, ",
                     "alcohol intake, physical activity, lipid-lowering medication use, glycemic status, ",
                     "hypertension status and total energy intake.")
note_fdr   <- paste0("FDR: Benjamini-Hochberg false discovery rate across the 8 UPF subclasses, ",
                     "calculated separately for each outcome.")
note_abbr  <- paste0("Abbreviations: CHD, coronary heart disease; CVD, cardiovascular disease; ",
                     "FDR, false discovery rate; MI, myocardial infarction; UPF, ultra-processed food.")
notes_S3   <- c(note_hr, note_main, note_fdr, note_abbr)

if (!all_same_cases) {
  notes_S3 <- append(notes_S3,
                     paste0("Cases / N are shown as a range because the number of participants with complete ",
                            "data differed slightly across subclass models."),
                     after = 1)
}

s3_num <- sub_raw %>%
  mutate(Outcome = factor(Outcome, levels = outcome_order)) %>%
  arrange(Outcome, Exposure_order) %>%
  select(Outcome, Exposure, Exposure_label, Exposure_order,
         N_used, N_events, beta, se, HR, Lower_CI, Upper_CI, P_value, q_value)

tabS3_long <- s3_num %>%
  transmute(key = unname(outcome_keys[as.character(Outcome)]),
            Exposure_order, Exposure_label,
            hr_ci = fmt_hr(HR, Lower_CI, Upper_CI),
            p     = fmt_p(P_value),
            fdr   = fmt_p(q_value))

tabS3_wide <- tabS3_long %>%
  pivot_wider(id_cols     = c(Exposure_order, Exposure_label),
              names_from  = key,
              values_from = c(hr_ci, p, fdr),
              names_glue  = "{key}_{.value}") %>%
  arrange(Exposure_order)

col_order <- c("Exposure_label",
               paste0(rep(unname(outcome_keys), each = 3), "_", c("hr_ci", "p", "fdr")))
tabS3 <- tabS3_wide %>% select(all_of(col_order))

header_top   <- c("", as.vector(rbind(unname(outcome_names[outcome_order]), "", "")))
header_cases <- c("", as.vector(rbind(cases_lab$label, "", "")))
header_sub   <- c("UPF subclass", rep(c("HR (95% CI)", "P", "FDR"), length(outcome_order)))

wb <- createWorkbook()
addWorksheet(wb, "Table S3")
n_col <- ncol(tabS3)

writeData(wb, "Table S3", t(header_top),   startRow = 1, startCol = 1, colNames = FALSE)
writeData(wb, "Table S3", t(header_cases), startRow = 2, startCol = 1, colNames = FALSE)
writeData(wb, "Table S3", t(header_sub),   startRow = 3, startCol = 1, colNames = FALSE)
writeData(wb, "Table S3", tabS3,           startRow = 4, startCol = 1, colNames = FALSE)

for (k in seq_along(outcome_order)) {
  cols_k <- (2 + 3 * (k - 1)):(4 + 3 * (k - 1))
  mergeCells(wb, "Table S3", cols = cols_k, rows = 1)
  mergeCells(wb, "Table S3", cols = cols_k, rows = 2)
}

body_end   <- 3 + nrow(tabS3)
note_start <- body_end + 2
writeData(wb, "Table S3", notes_S3, startRow = note_start, startCol = 1)

top_style   <- createStyle(textDecoration = "bold", halign = "center", valign = "center",
                           border = "bottom")
cases_style <- createStyle(textDecoration = "italic", halign = "center", valign = "center",
                           fontSize = 10)
sub_style   <- createStyle(textDecoration = "bold", halign = "center", valign = "center",
                           border = "bottom", wrapText = TRUE)
body_style  <- createStyle(halign = "center", valign = "center")
lab_style   <- createStyle(halign = "left",   valign = "center", wrapText = TRUE)
note_style  <- createStyle(wrapText = TRUE, valign = "top", fontSize = 10)
addStyle(wb, "Table S3", top_style,   rows = 1, cols = 2:n_col, gridExpand = TRUE)
addStyle(wb, "Table S3", cases_style, rows = 2, cols = 2:n_col, gridExpand = TRUE)
addStyle(wb, "Table S3", sub_style,   rows = 3, cols = 1:n_col, gridExpand = TRUE)
addStyle(wb, "Table S3", body_style,  rows = 4:body_end, cols = 2:n_col, gridExpand = TRUE)
addStyle(wb, "Table S3", lab_style,   rows = 4:body_end, cols = 1, gridExpand = TRUE)
addStyle(wb, "Table S3", createStyle(border = "bottom"), rows = body_end, cols = 1:n_col,
         gridExpand = TRUE, stack = TRUE)

note_rows <- note_start:(note_start + length(notes_S3) - 1)
for (r in note_rows) {
  mergeCells(wb, "Table S3", cols = 1:n_col, rows = r)
  addStyle(wb, "Table S3", note_style, rows = r, cols = 1:n_col, gridExpand = TRUE)
}
setRowHeights(wb, "Table S3", rows = note_rows, heights = 30)

setColWidths(wb, "Table S3", cols = 1:n_col,
             widths = c(34, rep(c(17, 7, 7), length(outcome_order))))

setwd(tabsFigs)
saveWorkbook(wb, "tableS3_UPFsubclass_cvdOutcomes.xlsx", overwrite = TRUE)
setwd(result)
write.csv(s3_num, "tableS3_UPFsubclass_cvdOutcomes_numbers.csv", row.names = FALSE)
