# Supplemental Tables S6 (baseline lipids) and S7 (5-year lipid change) with incident outcomes.
# Input: 3.1-3.3 and 4.1-4.3 results. Output: Tables S6 and S7 (result/tabsFigs).

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
raw_dir  <- "path/to/SHS_raw_data"
library(tidyverse)
library(openxlsx)
readyData <- file.path(proj_dir, "data")
result    <- file.path(proj_dir, "result")
tabsFigs  <- file.path(proj_dir, "result", "tabsFigs")
lipid_info_path <- file.path(raw_dir, "lipidomics data", "lipid_4_6_compare.csv")

if (!dir.exists(tabsFigs)) dir.create(tabsFigs, recursive = TRUE)

exposure_files <- c("longitu" = "Baseline", "change" = "Change")
outcome_files  <- c("cvd" = "CVD", "chd" = "CHD", "mi" = "MI", "stroke" = "Stroke")
set_files      <- c("riskLipids" = "Positive", "protLipids" = "Inverse")

setwd(result)
res_long <- expand_grid(ex  = names(exposure_files),
                        out = names(outcome_files),
                        st  = names(set_files)) %>%
  pmap_dfr(function(ex, out, st) {
    f <- paste0(ex, "_", out, "_", st, ".csv")
    d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    if (!"N_total"  %in% names(d)) d$N_total  <- NA
    if (!"N_events" %in% names(d)) d$N_events <- NA
    d %>%
      dplyr::select(lipids_name, beta, se, HR, lower_CI, upper_CI, pvals, qvals,
                    N_total_file = N_total, N_events_file = N_events) %>%
      mutate(exposure  = exposure_files[[ex]],
             outcome   = outcome_files[[out]],
             Direction = set_files[[st]],
             file      = f)
  })

key_lipids_risk <- read.table("UPF_lipids_risk.csv", sep = ",", header = TRUE, check.names = FALSE)
key_lipids_prot <- read.table("UPF_lipids_prot.csv", sep = ",", header = TRUE, check.names = FALSE)

f_fig <- "lipids42_cvdOutcomes_bsl_chg_long.csv"
if (file.exists(f_fig)) {
  fig_long <- read.csv(f_fig, stringsAsFactors = FALSE, check.names = FALSE)
  chk_fig  <- res_long %>%
    dplyr::select(lipids_name, exposure, outcome, HR, pvals) %>%
    full_join(fig_long %>% dplyr::select(lipids_name, exposure, outcome, HR, pvals),
              by = c("lipids_name", "exposure", "outcome"), suffix = c("_s6", "_fig"))
} else {
}

lipid_info <- read.csv(lipid_info_path, header = TRUE, stringsAsFactors = FALSE,
                       check.names = FALSE) %>%
  dplyr::select(lipids_name, Class_short) %>%
  distinct(lipids_name, .keep_all = TRUE)

lip42 <- res_long %>%
  distinct(Direction, lipids_name) %>%
  left_join(lipid_info, by = "lipids_name") %>%
  mutate(Class_short = ifelse(is.na(Class_short) | Class_short == "", "N/A", Class_short))
lip42[lip42$Class_short == "N/A",                       "Class_short"] <- "Unknown"
lip42[lip42$Class_short %in% c("FAHFA", "scFA", "USFA"), "Class_short"] <- "FA"
lip42[lip42$Class_short %in% c("Glc"),                  "Class_short"] <- "GlcCer"

cls_order <- c("AC", "CE", "CER", "DAG", "FA", "GlcCer", "PC", "PE", "PG", "PI",
               "PS", "SM", "TAG", "Unknown")

lip42 <- lip42 %>%
  mutate(Direction   = factor(Direction, levels = c("Positive", "Inverse")),
         Class_short = factor(Class_short, levels = cls_order),
         name_rank   = order(str_order(lipids_name, numeric = TRUE))) %>%
  arrange(Direction, Class_short, name_rank) %>%
  mutate(row_order = row_number()) %>%
  dplyr::select(Direction, Class_short, lipids_name, row_order)

setwd(readyData)
load("lipidomics_p4.Rda")
load("dat_AimOne_p4.Rda")
load("outcome_mi_stroke_p4.Rda")
load("dat_lipidChange_P5landmark.Rda")

dt_bsl    <- as.data.frame(lipidomics_p4) %>%
  inner_join(dat_AimOne_p4, by = c("ID" = "IDNO"))
dt_bsl_ms <- dt_bsl %>%
  inner_join(outcome_mi_stroke_p4, by = c("ID" = "IDNO"))

covs_bsl <- c("S4AGE", "SEX", "CENTER", "S4EDU", "S4BMI", "S4SMOKE", "S4ETOH", "PDM3",
              "LIPIDLOW", "S4ADADM2020", "S4USHTN", "S4DT_KCAL", "FID")
covs_chg <- c("S5AGE", "SEX", "CENTER", "S5EDU", "S5BMI", "S5SMOKE", "S5ETOH", "ped_avg3day",
              "LIPIDLOW_p5", "S5ADADM2020", "S5USHTN", "S5DT_KCAL", "BMI_diff", "DT_KCAL_diff", "FID")

samp_def <- list(
  list(exposure = "Baseline", outcome = "CVD",    d = dt_bsl,    time = "TIME_cvd",     event = "cvd",             atrisk = NA),
  list(exposure = "Baseline", outcome = "CHD",    d = dt_bsl,    time = "TIME_chd",     event = "chd",             atrisk = NA),
  list(exposure = "Baseline", outcome = "MI",     d = dt_bsl_ms, time = "TIME_mi",      event = "mi",              atrisk = NA),
  list(exposure = "Baseline", outcome = "Stroke", d = dt_bsl_ms, time = "TIME_stroke",  event = "definite_stroke", atrisk = NA),
  list(exposure = "Change",   outcome = "CVD",    d = dat_lipidChange, time = "TIME_cvd5",    event = "cvd5",    atrisk = "atrisk_cvd5"),
  list(exposure = "Change",   outcome = "CHD",    d = dat_lipidChange, time = "TIME_chd5",    event = "chd5",    atrisk = "atrisk_chd5"),
  list(exposure = "Change",   outcome = "MI",     d = dat_lipidChange, time = "TIME_mi5",     event = "mi5",     atrisk = "atrisk_mi5"),
  list(exposure = "Change",   outcome = "Stroke", d = dat_lipidChange, time = "TIME_stroke5", event = "stroke5", atrisk = "atrisk_stroke5"))

n_event <- function(x) sum(x %in% c(1, TRUE, "1", "TRUE"))

n_long <- map_dfr(samp_def, function(s) {
  d <- as.data.frame(s$d)
  if (!is.na(s$atrisk)) d <- d[d[[s$atrisk]] %in% TRUE, ]
  map_dfr(lip42$lipids_name, function(lp) {
    lip_vars <- if (s$exposure == "Baseline") lp else paste0(lp, c("_chg", "_bsl"))
    vars     <- c(lip_vars, s$time, s$event, if (s$exposure == "Baseline") covs_bsl else covs_chg)
    miss     <- setdiff(vars, names(d))
    cc <- complete.cases(d[, vars])
    tibble(exposure = s$exposure, outcome = s$outcome, lipids_name = lp,
           N = sum(cc), N_cases = n_event(d[[s$event]][cc]))
  })
})

n_model <- n_long %>%
  group_by(exposure, outcome) %>%
  summarise(N_min = min(N), N_max = max(N), cases_min = min(N_cases), cases_max = max(N_cases),
            .groups = "drop")

fmt_hr  <- function(h, lo, hi) sprintf("%.2f (%.2f, %.2f)", h, lo, hi)
fmt_sci <- function(p) formatC(p, format = "E", digits = 1)
fmt_n   <- function(x) format(x, big.mark = ",", trim = TRUE)
fmt_rng <- function(a, b) ifelse(a == b, fmt_n(a), paste0(fmt_n(a), "\u2013", fmt_n(b)))

n_model <- n_model %>%
  mutate(head_lab = paste0(outcome, "\n(cases/N = ", fmt_rng(cases_min, cases_max), "/",
                           fmt_rng(N_min, N_max), ")"))
n_varies <- any(n_model$N_min != n_model$N_max | n_model$cases_min != n_model$cases_max)

cls_full <- c(AC = "acylcarnitine", CE = "cholesteryl ester", CER = "ceramide",
              DAG = "diacylglycerol", FA = "fatty acid", GlcCer = "glucosylceramide",
              PC = "phosphatidylcholine", PE = "phosphatidylethanolamine",
              PG = "phosphatidylglycerol", PI = "phosphatidylinositol",
              PS = "phosphatidylserine", SM = "sphingomyelin", TAG = "triacylglycerol")
cls_in   <- intersect(names(cls_full), as.character(unique(lip42$Class_short)))
abbr_all <- c(paste0(cls_in, ", ", cls_full[cls_in]),
              "CHD, coronary heart disease", "CI, confidence interval",
              "CVD, cardiovascular disease", "FDR, false discovery rate", "HR, hazard ratio",
              "MI, myocardial infarction", "UPF, ultra-processed food")
abbr_all <- abbr_all[order(toupper(abbr_all))]
note_abbr <- paste0("Abbreviations: ", paste(abbr_all, collapse = "; "), ".")

note_fdr   <- paste0("FDR q: Benjamini-Hochberg false discovery rate across the 42 UPF-related lipids, ",
                     "within each outcome. Bold: P < 0.05.")
note_set   <- paste0("The 42 lipids are the named lipids associated with UPF intake at FDR < 0.05 ",
                     "(Supplemental Table S4). Positive markers (higher with higher UPF intake, n = ",
                     sum(lip42$Direction == "Positive"), ") are listed first, then inverse markers (n = ",
                     sum(lip42$Direction == "Inverse"), "); within each group, lipids are ordered by ",
                     "class and then by name.")
note_out   <- paste0("CVD, fatal and nonfatal cardiovascular disease; CHD, coronary heart disease; ",
                     "MI, myocardial infarction; Stroke, definite stroke. Cases/N: number of incident ",
                     "cases / number of participants in the model (complete cases)",
                     ifelse(n_varies, "; a range is shown where the number differs across lipids.", "."))

note_hr_a  <- paste0("Hazard ratios (HRs) and 95% CIs are per 1-SD higher baseline (2001-2003) lipid ",
                     "level, from Cox proportional hazards models with a family-level frailty term. ",
                     "Follow-up ran from the baseline examination to the event, death or ",
                     "December 31, 2020, whichever came first.")
note_adj_a <- paste0("Models were adjusted for baseline age, sex, study center, education, body mass ",
                     "index, smoking, alcohol intake, physical activity, lipid-lowering medication use, ",
                     "glycemic status, hypertension and total energy intake.")
notes_a    <- c(note_hr_a, note_adj_a, note_fdr, note_set, note_out, note_abbr)

note_hr_b  <- paste0("HRs and 95% CIs are per 1-SD higher 5-year change in lipid level (follow-up ",
                     "[2006-2009] minus baseline [2001-2003]), from landmark Cox proportional hazards ",
                     "models with a family-level frailty term. Follow-up started at the follow-up ",
                     "examination and included only participants free of the given outcome at that ",
                     "examination; it ran to the event, death or December 31, 2020, whichever came first.")
note_adj_b <- paste0("Models were adjusted for age, sex, study center, education, body mass index, ",
                     "smoking, alcohol intake, physical activity, lipid-lowering medication use, ",
                     "glycemic status, hypertension and total energy intake at the follow-up ",
                     "examination, plus 5-year changes in body mass index and total energy intake and ",
                     "the baseline level of the lipid. Associations of baseline lipid levels are ",
                     "shown in Supplemental Table S6.")
notes_b    <- c(note_hr_b, note_adj_b, note_fdr, note_set, note_out, note_abbr)

title_a <- paste0("Supplemental Table S6. Associations of baseline levels of the 42 UPF-related ",
                  "lipids with incident CVD, CHD, MI and stroke")
title_b <- paste0("Supplemental Table S7. Associations of 5-year changes in the 42 UPF-related ",
                  "lipids with incident CVD, CHD, MI and stroke")

s6_num <- res_long %>%
  dplyr::select(-Direction) %>%
  inner_join(lip42, by = "lipids_name") %>%
  left_join(n_long, by = c("exposure", "outcome", "lipids_name")) %>%
  mutate(exposure = factor(exposure, levels = c("Baseline", "Change")),
         outcome  = factor(outcome,  levels = outcome_files)) %>%
  arrange(exposure, row_order, outcome) %>%
  dplyr::select(exposure, Direction, Class_short, lipids_name, outcome, N, N_cases,
                beta, se, HR, lower_CI, upper_CI, pvals, qvals, file)

make_wide <- function(ex) {
  d <- s6_num %>% filter(exposure == ex)
  tab <- d %>%
    transmute(Direction, row_order = match(lipids_name, lip42$lipids_name),
              `Lipid class` = as.character(Class_short), Lipid = lipids_name, outcome,
              hr = fmt_hr(HR, lower_CI, upper_CI), p = fmt_sci(pvals), q = fmt_sci(qvals)) %>%
    pivot_wider(names_from = outcome, values_from = c(hr, p, q),
                names_glue = "{outcome}_{.value}") %>%
    arrange(row_order)
  col_out <- unname(unlist(lapply(outcome_files, function(o) paste0(o, c("_hr", "_p", "_q")))))
  sig <- d %>%
    transmute(lipids_name, outcome, sig = pvals < 0.05) %>%
    pivot_wider(names_from = outcome, values_from = sig)
  sig <- sig[match(tab$Lipid, sig$lipids_name), unname(outcome_files)]
  list(tab = tab %>% dplyr::select(Direction, `Lipid class`, Lipid, all_of(col_out)),
       sig = as.matrix(sig))
}
wide_a <- make_wide("Baseline")
wide_b <- make_wide("Change")

title_style  <- createStyle(textDecoration = "bold", halign = "left", valign = "center")
header_style <- createStyle(textDecoration = "bold", halign = "center", valign = "center",
                            wrapText = TRUE)
sub_style    <- createStyle(textDecoration = c("bold", "italic"), halign = "left", valign = "center")
left_style   <- createStyle(halign = "left",   valign = "center")
num_style    <- createStyle(halign = "center", valign = "center")
bold_style   <- createStyle(textDecoration = "bold")
note_style   <- createStyle(wrapText = TRUE, valign = "top", fontSize = 10)

write_s6_sheet <- function(wb, sh, w, title, notes, ex) {
  addWorksheet(wb, sh)
  tab   <- w$tab
  n_col <- ncol(tab) - 1
  n_out <- length(outcome_files)
  n_pos         <- sum(tab$Direction == "Positive")
  n_inv         <- sum(tab$Direction == "Inverse")
  row_sub_pos   <- 4
  row_pos_start <- 5
  row_pos_end   <- row_pos_start + n_pos - 1
  row_sub_inv   <- row_pos_end + 1
  row_inv_start <- row_sub_inv + 1
  row_inv_end   <- row_inv_start + n_inv - 1
  note_start    <- row_inv_end + 2
  writeData(wb, sh, title, startRow = 1, startCol = 1)
  mergeCells(wb, sh, cols = 1:n_col, rows = 1)
  addStyle(wb, sh, title_style, rows = 1, cols = 1:n_col, gridExpand = TRUE)
  head1 <- n_model %>% filter(exposure == ex) %>%
    mutate(outcome = factor(outcome, levels = outcome_files)) %>% arrange(outcome)
  writeData(wb, sh, t(c("Lipid class", "Lipid")), startRow = 2, startCol = 1, colNames = FALSE)
  for (k in 1:n_out) {
    c1 <- 2 + 3 * (k - 1) + 1
    writeData(wb, sh, head1$head_lab[k], startRow = 2, startCol = c1)
    mergeCells(wb, sh, cols = c1:(c1 + 2), rows = 2)
  }
  mergeCells(wb, sh, cols = 1, rows = 2:3)
  mergeCells(wb, sh, cols = 2, rows = 2:3)
  writeData(wb, sh, t(rep(c("HR (95% CI)", "P", "FDR q"), n_out)), startRow = 3, startCol = 3,
            colNames = FALSE)
  addStyle(wb, sh, header_style, rows = 2:3, cols = 1:n_col, gridExpand = TRUE)
  addStyle(wb, sh, createStyle(border = "top"),    rows = 2, cols = 1:n_col, gridExpand = TRUE, stack = TRUE)
  addStyle(wb, sh, createStyle(border = "bottom"), rows = 2, cols = 3:n_col, gridExpand = TRUE, stack = TRUE)
  addStyle(wb, sh, createStyle(border = "bottom"), rows = 3, cols = 1:n_col, gridExpand = TRUE, stack = TRUE)
  setRowHeights(wb, sh, rows = 2, heights = 30)
  writeData(wb, sh, paste0("Positively associated with UPF intake (n = ", n_pos, ")"),
            startRow = row_sub_pos, startCol = 1)
  writeData(wb, sh, tab %>% filter(Direction == "Positive") %>% dplyr::select(-Direction),
            startRow = row_pos_start, startCol = 1, colNames = FALSE)
  writeData(wb, sh, paste0("Inversely associated with UPF intake (n = ", n_inv, ")"),
            startRow = row_sub_inv, startCol = 1)
  writeData(wb, sh, tab %>% filter(Direction == "Inverse") %>% dplyr::select(-Direction),
            startRow = row_inv_start, startCol = 1, colNames = FALSE)
  for (r in c(row_sub_pos, row_sub_inv)) {
    mergeCells(wb, sh, cols = 1:n_col, rows = r)
    addStyle(wb, sh, sub_style, rows = r, cols = 1:n_col, gridExpand = TRUE)
  }
  body_rows <- c(row_pos_start:row_pos_end, row_inv_start:row_inv_end)
  addStyle(wb, sh, left_style, rows = body_rows, cols = 1:2,     gridExpand = TRUE)
  addStyle(wb, sh, num_style,  rows = body_rows, cols = 3:n_col, gridExpand = TRUE)
  for (i in seq_along(body_rows)) {
    for (k in 1:n_out) {
      if (w$sig[i, k]) {
        c1 <- 2 + 3 * (k - 1) + 1
        addStyle(wb, sh, bold_style, rows = body_rows[i], cols = c1:(c1 + 2),
                 gridExpand = TRUE, stack = TRUE)
      }
    }
  }
  addStyle(wb, sh, createStyle(border = "bottom"), rows = row_inv_end, cols = 1:n_col,
           gridExpand = TRUE, stack = TRUE)
  writeData(wb, sh, notes, startRow = note_start, startCol = 1)
  for (r in note_start:(note_start + length(notes) - 1)) {
    mergeCells(wb, sh, cols = 1:n_col, rows = r)
    addStyle(wb, sh, note_style, rows = r, cols = 1:n_col, gridExpand = TRUE)
  }
  setRowHeights(wb, sh, rows = note_start:(note_start + length(notes) - 1), heights = 30)
  setColWidths(wb, sh, cols = 1:n_col, widths = c(10, 26, rep(c(17, 9, 9), n_out)))
  freezePane(wb, sh, firstActiveRow = 4, firstActiveCol = 3)
  pageSetup(wb, sh, orientation = "landscape", fitToWidth = TRUE, fitToHeight = FALSE)
  invisible(wb)
}

wb_s6 <- createWorkbook()
write_s6_sheet(wb_s6, "Table S6", wide_a, title_a, notes_a, "Baseline")
wb_s7 <- createWorkbook()
write_s6_sheet(wb_s7, "Table S7", wide_b, title_b, notes_b, "Change")

setwd(tabsFigs)
saveWorkbook(wb_s6, "tableS6_lipids42_cvdOutcomes_baseline.xlsx", overwrite = TRUE)
saveWorkbook(wb_s7, "tableS7_lipids42_cvdOutcomes_change.xlsx",   overwrite = TRUE)
setwd(result)
write.csv(s6_num, "tableS6S7_lipids42_cvdOutcomes_numbers.csv", row.names = FALSE)
