# Supplemental Table S8: mediation by 5-year lipid change, main and sensitivity models.
# Input: 5.6 results. Output: Table S8 (result/tabsFigs).

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

n_sims <- 1000
p_bold <- 0.05

setwd(result)
model_files   <- c(main           = "mediation_change_mediationPkg_survreg_main.csv",
                   sensNoBMIDMHTN = "mediation_change_mediationPkg_survreg_sensNoBMIDMHTN.csv")
model_labels  <- c(main           = "A. Main model",
                   sensNoBMIDMHTN = paste0("B. Sensitivity model: without adjustment for BMI, ",
                                           "glycemic status and hypertension"))
outcome_order <- c("Incident CVD", "Incident CHD", "Incident MI")

med_long <- imap_dfr(model_files, function(f, m) {
  read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) %>%
    mutate(model = m, file = f)
})

f_fig <- "fig_lipidsCVD_panel_b_mediationPkg_survreg_main_plotdata.csv"
if (file.exists(f_fig)) {
  fig_b   <- read.csv(f_fig, stringsAsFactors = FALSE, check.names = FALSE)
  chk_fig <- med_long %>%
    filter(model == "main", ACME_p < 0.05) %>%
    dplyr::select(outcome, lipids_name, PropMed, ACME_p) %>%
    full_join(fig_b %>% dplyr::select(outcome, lipids_name, PropMed, ACME_p),
              by = c("outcome", "lipids_name"), suffix = c("_s8", "_fig"))
} else {
}

n_p0 <- sum(med_long$ACME_p == 0)

lipid_info <- read.csv(lipid_info_path, header = TRUE, stringsAsFactors = FALSE,
                       check.names = FALSE) %>%
  dplyr::select(lipids_name, Class_short) %>%
  distinct(lipids_name, .keep_all = TRUE)

lip_med <- med_long %>%
  distinct(lipids_name) %>%
  left_join(lipid_info, by = "lipids_name") %>%
  mutate(Class_short = ifelse(is.na(Class_short) | Class_short == "", "N/A", Class_short))
lip_med[lip_med$Class_short == "N/A",                       "Class_short"] <- "Unknown"
lip_med[lip_med$Class_short %in% c("FAHFA", "scFA", "USFA"), "Class_short"] <- "FA"
lip_med[lip_med$Class_short %in% c("Glc"),                  "Class_short"] <- "GlcCer"

cls_order <- c("AC", "CE", "CER", "DAG", "FA", "GlcCer", "PC", "PE", "PG", "PI",
               "PS", "SM", "TAG", "Unknown")

lip_med <- lip_med %>%
  mutate(Class_short = factor(Class_short, levels = cls_order),
         name_rank   = order(str_order(lipids_name, numeric = TRUE)))

fmt_num <- function(x, d) {
  x <- round(x, d); x[x == 0] <- 0
  formatC(x, format = "f", digits = d)
}
fmt_ci  <- function(est, lo, hi, d) paste0(fmt_num(est, d), " (", fmt_num(lo, d), ", ", fmt_num(hi, d), ")")
fmt_pct <- function(x) ifelse(abs(x) < 0.1 & round(x, 2) != 0, fmt_num(x, 2), fmt_num(x, 1))
p_min   <- 2 / n_sims
fmt_sci <- function(p) ifelse(p == 0, paste0("<", formatC(p_min, format = "E", digits = 1)),
                              formatC(p, format = "E", digits = 1))
fmt_n   <- function(x) format(x, big.mark = ",", trim = TRUE)
fmt_rng <- function(a, b) ifelse(a == b, fmt_n(a), paste0(fmt_n(a), "\u2013", fmt_n(b)))
fmt_rng_y <- function(a, b) ifelse(fmt_num(a, 1) == fmt_num(b, 1), fmt_num(a, 1),
                                   paste0(fmt_num(a, 1), "\u2013", fmt_num(b, 1)))

cls_full <- c(AC = "acylcarnitine", CE = "cholesteryl ester", CER = "ceramide",
              DAG = "diacylglycerol", FA = "fatty acid", GlcCer = "glucosylceramide",
              PC = "phosphatidylcholine", PE = "phosphatidylethanolamine",
              PG = "phosphatidylglycerol", PI = "phosphatidylinositol",
              PS = "phosphatidylserine", SM = "sphingomyelin", TAG = "triacylglycerol")
cls_in   <- intersect(names(cls_full), as.character(unique(lip_med$Class_short)))
abbr_all <- c(paste0(cls_in, ", ", cls_full[cls_in]),
              "BMI, body mass index", "CHD, coronary heart disease", "CI, confidence interval",
              "CVD, cardiovascular disease", "FDR, false discovery rate",
              "MI, myocardial infarction", "TR, time ratio", "UPF, ultra-processed food")
abbr_all <- abbr_all[order(toupper(abbr_all))]
note_abbr <- paste0("Abbreviations: ", paste(abbr_all, collapse = "; "), ".")

note_design <- paste0("Mediation analyses used the landmark sample: follow-up started at the follow-up ",
                      "examination (2006-2009) and included only participants free of the given outcome ",
                      "at that examination; it ran to the event, death or December 31, 2020, whichever ",
                      "came first. The exposure was baseline (2001-2003) UPF intake (% of energy) and the ",
                      "mediator was the 5-year change in lipid level (follow-up minus baseline, per SD).")
note_cand   <- paste0("Candidate mediators for each outcome were the UPF-related lipids whose 5-year ",
                      "change was associated with that outcome at P < 0.05 (Supplemental Table S7); all ",
                      "candidates tested are shown. Stroke was not analyzed because UPF intake was not ",
                      "associated with incident stroke.")
note_paths  <- paste0("a path: \u03b2 (SD of 5-year lipid change per 10% higher energy from UPF) from ",
                      "linear regression. b and c\u2032 paths: TRs per 1-SD higher 5-year lipid change (b) ",
                      "and per 10% higher energy from UPF (c\u2032, direct effect) from Weibull accelerated ",
                      "failure time models that included both UPF intake and lipid change; TR < 1 means a ",
                      "shorter time to the event (higher risk). 95% CIs use family-clustered robust ",
                      "standard errors.")
note_ind    <- paste0("Indirect effect: average causal mediation effect from the R mediation package ",
                      "(quasi-Bayesian Monte Carlo approximation, ", fmt_n(n_sims), " simulations), ",
                      "comparing the median UPF intake with the median + 10% of energy. % mediated = ",
                      "indirect effect / total effect, with simulation-based 95% CI; wide CIs reflect ",
                      "total effects close to the null in the landmark sample. Inconsistent mediation: ",
                      "the indirect and direct effects had opposite signs, so % mediated is not ",
                      "interpretable.",
                      ifelse(n_p0 > 0, paste0(" P values below the simulation resolution are shown as < ",
                                              formatC(p_min, format = "E", digits = 1), "."), ""))
note_fdr    <- paste0("FDR q: Benjamini-Hochberg false discovery rate across the candidate mediators ",
                      "within each outcome and model. Bold: indirect effect P < ", p_bold, ".")
note_adj    <- paste0("Main model: adjusted for baseline age, sex, study center, education, BMI, smoking, ",
                      "alcohol intake, physical activity, lipid-lowering medication use, glycemic status, ",
                      "hypertension and total energy intake, and the baseline level of the lipid (in both ",
                      "the mediator and outcome models). Sensitivity model: same, without BMI, glycemic ",
                      "status and hypertension, which may lie on the pathway from UPF intake to CVD; both ",
                      "models included the same participants.")
note_out    <- paste0("CVD, fatal and nonfatal cardiovascular disease; CHD, coronary heart disease; MI, ",
                      "myocardial infarction. Cases/N: number of incident cases / number of participants ",
                      "in the model (complete cases); a range is shown where the number differs across ",
                      "lipids. Median follow-up is from the follow-up examination.")
notes_s8    <- c(note_design, note_cand, note_adj, note_paths, note_ind, note_fdr, note_out, note_abbr)

title_s8 <- paste0("Supplemental Table S8. Mediation of the associations of UPF intake with incident ",
                   "CVD, CHD and MI by 5-year changes in UPF-related lipids")

tab_num <- med_long %>%
  inner_join(lip_med, by = "lipids_name") %>%
  mutate(model   = factor(model,   levels = names(model_files)),
         outcome = factor(outcome, levels = outcome_order)) %>%
  arrange(model, outcome, Class_short, name_rank)

tab_fmt <- tab_num %>%
  mutate(`Lipid class` = as.character(Class_short),
         Lipid  = lipids_name,
         a_txt  = fmt_ci(a_beta, a_lo, a_hi, 3),
         b_txt  = fmt_ci(b_TR, b_TR_lo, b_TR_hi, 2),
         c_txt  = fmt_ci(c_TR, c_TR_lo, c_TR_hi, 2),
         pm_txt = ifelse(same_sign,
                         paste0(fmt_pct(100 * PropMed), " (", fmt_pct(100 * PropMed_lo), ", ",
                                fmt_pct(100 * PropMed_hi), ")"),
                         "Inconsistent mediation"),
         p_txt  = fmt_sci(ACME_p),
         q_txt  = fmt_sci(ACME_q),
         sig    = ACME_p < p_bold)

sub_info <- tab_num %>%
  group_by(model, outcome) %>%
  summarise(n_med  = n(),
            N_min  = min(N_total),  N_max  = max(N_total),
            ev_min = min(N_events), ev_max = max(N_events),
            fu_min = min(fu_years_median), fu_max = max(fu_years_median),
            .groups = "drop") %>%
  mutate(sub_lab = paste0(outcome, " (cases/N = ", fmt_rng(ev_min, ev_max), "/", fmt_rng(N_min, N_max),
                          "; median follow-up ", fmt_rng_y(fu_min, fu_max), " years; ",
                          n_med, " candidate mediators)"))

lay <- map_dfr(names(model_files), function(m) {
  blk  <- tibble(type = "block", `Lipid class` = model_labels[[m]])
  subs <- map_dfr(outcome_order, function(o) {
    d <- tab_fmt %>% filter(model == m, outcome == o)
    bind_rows(tibble(type = "sub",
                     `Lipid class` = sub_info$sub_lab[sub_info$model == m & sub_info$outcome == o]),
              d %>% transmute(type = "data", `Lipid class`, Lipid, a_txt, b_txt, c_txt,
                              pm_txt, p_txt, q_txt, sig))
  })
  bind_rows(blk, subs)
})

title_style  <- createStyle(textDecoration = "bold", halign = "left", valign = "center")
header_style <- createStyle(textDecoration = "bold", halign = "center", valign = "center",
                            wrapText = TRUE)
block_style  <- createStyle(textDecoration = "bold", halign = "left", valign = "center",
                            fgFill = "#F2F2F2")
sub_style    <- createStyle(textDecoration = c("bold", "italic"), halign = "left", valign = "center")
left_style   <- createStyle(halign = "left",   valign = "center")
num_style    <- createStyle(halign = "center", valign = "center")
bold_style   <- createStyle(textDecoration = "bold")
note_style   <- createStyle(wrapText = TRUE, valign = "top", fontSize = 10)

wb    <- createWorkbook()
sh    <- "Table S8"
addWorksheet(wb, sh)
n_col      <- 8
body_start <- 4
body_rows  <- body_start + seq_len(nrow(lay)) - 1
note_start <- max(body_rows) + 2

writeData(wb, sh, title_s8, startRow = 1, startCol = 1)
mergeCells(wb, sh, cols = 1:n_col, rows = 1)
addStyle(wb, sh, title_style, rows = 1, cols = 1:n_col, gridExpand = TRUE)

writeData(wb, sh, t(c("Lipid class", "Lipid",
                      "a path\nUPF \u2192 \u0394lipid",
                      "b path\n\u0394lipid \u2192 outcome",
                      "c\u2032 path (direct effect)\nUPF \u2192 outcome",
                      "Indirect effect\nUPF \u2192 \u0394lipid \u2192 outcome", "", "")),
          startRow = 2, startCol = 1, colNames = FALSE)
mergeCells(wb, sh, cols = 1,   rows = 2:3)
mergeCells(wb, sh, cols = 2,   rows = 2:3)
mergeCells(wb, sh, cols = 6:8, rows = 2)
writeData(wb, sh, t(c("\u03b2 (95% CI)", "TR (95% CI)", "TR (95% CI)",
                      "% mediated (95% CI)", "P", "FDR q")),
          startRow = 3, startCol = 3, colNames = FALSE)
addStyle(wb, sh, header_style, rows = 2:3, cols = 1:n_col, gridExpand = TRUE)
addStyle(wb, sh, createStyle(border = "top"),    rows = 2, cols = 1:n_col, gridExpand = TRUE, stack = TRUE)
addStyle(wb, sh, createStyle(border = "bottom"), rows = 2, cols = 6:8,     gridExpand = TRUE, stack = TRUE)
addStyle(wb, sh, createStyle(border = "bottom"), rows = 3, cols = 1:n_col, gridExpand = TRUE, stack = TRUE)
setRowHeights(wb, sh, rows = 2, heights = 30)

writeData(wb, sh, lay %>% dplyr::select(`Lipid class`, Lipid, a_txt, b_txt, c_txt, pm_txt, p_txt, q_txt),
          startRow = body_start, startCol = 1, colNames = FALSE)

row_block <- body_rows[lay$type == "block"]
row_sub   <- body_rows[lay$type == "sub"]
row_data  <- body_rows[lay$type == "data"]
for (r in row_block) {
  mergeCells(wb, sh, cols = 1:n_col, rows = r)
  addStyle(wb, sh, block_style, rows = r, cols = 1:n_col, gridExpand = TRUE)
}
for (r in row_sub) {
  mergeCells(wb, sh, cols = 1:n_col, rows = r)
  addStyle(wb, sh, sub_style, rows = r, cols = 1:n_col, gridExpand = TRUE)
}
addStyle(wb, sh, left_style, rows = row_data, cols = 1:2,     gridExpand = TRUE)
addStyle(wb, sh, num_style,  rows = row_data, cols = 3:n_col, gridExpand = TRUE)
addStyle(wb, sh, createStyle(border = "top"), rows = row_block[-1], cols = 1:n_col,
         gridExpand = TRUE, stack = TRUE)

row_sig <- body_rows[lay$sig %in% TRUE]
if (length(row_sig) > 0) {
  addStyle(wb, sh, bold_style, rows = row_sig, cols = 6:8, gridExpand = TRUE, stack = TRUE)
}

addStyle(wb, sh, createStyle(border = "bottom"), rows = max(body_rows), cols = 1:n_col,
         gridExpand = TRUE, stack = TRUE)

writeData(wb, sh, notes_s8, startRow = note_start, startCol = 1)
for (r in note_start:(note_start + length(notes_s8) - 1)) {
  mergeCells(wb, sh, cols = 1:n_col, rows = r)
  addStyle(wb, sh, note_style, rows = r, cols = 1:n_col, gridExpand = TRUE)
}
setRowHeights(wb, sh, rows = note_start:(note_start + length(notes_s8) - 1), heights = 42)

setColWidths(wb, sh, cols = 1:n_col, widths = c(10, 22, 24, 18, 22, 26, 10, 10))
freezePane(wb, sh, firstActiveRow = body_start, firstActiveCol = 3)
pageSetup(wb, sh, orientation = "landscape", fitToWidth = TRUE, fitToHeight = FALSE)

setwd(tabsFigs)
saveWorkbook(wb, "tableS8_mediation_lipidChange.xlsx", overwrite = TRUE)
setwd(result)
write.csv(tab_num %>%
            dplyr::select(model, outcome, Class_short, lipids_name, N_total, N_events, fu_years_median,
                          upf_contrast, a_beta, a_lo, a_hi, a_p, b_TR, b_TR_lo, b_TR_hi, b_p,
                          c_TR, c_TR_lo, c_TR_hi, ACME_p, ACME_q,
                          PropMed, PropMed_lo, PropMed_hi, same_sign, covariates, file),
          "tableS8_mediation_lipidChange_numbers.csv", row.names = FALSE)
