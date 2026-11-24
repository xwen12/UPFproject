# Supplemental Table S5: center-split sensitivity analysis (discovery AZ + SD, replication OK).
# Input: result/sensCenter_AZSD_OK_combined_all.csv (2.4). Output: Table S5 (result/tabsFigs).

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
raw_dir  <- "path/to/SHS_raw_data"
library(tidyverse)
library(openxlsx)
result   <- file.path(proj_dir, "result")
tabsFigs <- file.path(proj_dir, "result", "tabsFigs")
lipid_info_path <- file.path(raw_dir, "lipidomics data", "lipid_4_6_compare.csv")

if (!dir.exists(tabsFigs)) dir.create(tabsFigs, recursive = TRUE)

setwd(result)
comb <- read.table("sensCenter_AZSD_OK_combined_all.csv", sep = ",", header = TRUE, check.names = FALSE)
comb <- comb %>%
  mutate(across(c(known, in42, disc_hit, repl_primary), as.logical))

upf_lipids42 <- read.table("UPF_lipids_all.csv", sep = ",", header = TRUE, check.names = FALSE)

s5_raw <- comb %>%
  filter(disc_hit, known)
n_known_tested <- sum(comb$known & !is.na(comb$p_disc))

n_hits     <- nrow(s5_raw)
n_het05    <- sum(s5_raw$p_het < 0.05)
med_I2     <- median(s5_raw$I2_het)

lipid_info <- read.csv(lipid_info_path, header = TRUE, stringsAsFactors = FALSE,
                       check.names = FALSE)

lipid_info <- lipid_info %>%
  select(lipids_name, Class_short) %>%
  distinct(lipids_name, .keep_all = TRUE)

s5_dt <- s5_raw %>%
  left_join(lipid_info, by = "lipids_name")

s5_dt <- s5_dt %>%
  mutate(Class_short = ifelse(is.na(Class_short) | Class_short == "", "N/A", Class_short))
s5_dt[s5_dt$Class_short == "N/A",                       "Class_short"] <- "Unknown"
s5_dt[s5_dt$Class_short %in% c("FAHFA", "scFA", "USFA"), "Class_short"] <- "FA"
s5_dt[s5_dt$Class_short %in% c("Glc"),                  "Class_short"] <- "GlcCer"

cls_order <- c("AC", "CE", "CER", "DAG", "FA", "GlcCer", "PC", "PE", "PG", "PI",
               "PS", "SM", "TAG", "Unknown")
cls_extra <- sort(setdiff(unique(s5_dt$Class_short), cls_order))
cls_order <- c(setdiff(cls_order, "Unknown"), cls_extra, "Unknown")

s5_num <- s5_dt %>%
  mutate(Direction     = factor(ifelse(beta_disc > 0, "Positive", "Inverse"),
                                levels = c("Positive", "Inverse")),
         Class_short   = factor(Class_short, levels = cls_order),
         name_rank     = order(str_order(lipids_name, numeric = TRUE)),
         lower_CI_disc = beta_disc - 1.96 * se_disc,
         upper_CI_disc = beta_disc + 1.96 * se_disc,
         lower_CI_repl = beta_repl - 1.96 * se_repl,
         upper_CI_repl = beta_repl + 1.96 * se_repl,
         lower_CI_meta = beta_meta - 1.96 * se_meta,
         upper_CI_meta = beta_meta + 1.96 * se_meta) %>%
  arrange(Direction, Class_short, name_rank) %>%
  select(Direction, Class_short, lipids_name,
         beta_disc, se_disc, lower_CI_disc, upper_CI_disc, p_disc, N_disc,
         beta_repl, se_repl, lower_CI_repl, upper_CI_repl, p_repl, N_repl,
         beta_meta, se_meta, lower_CI_meta, upper_CI_meta, p_meta,
         Q_het, I2_het, p_het,
         replicated = repl_primary, in42, beta_full, qvals_full)

fmt_beta <- function(b, lo, hi) sprintf("%.3f (%.3f, %.3f)", round(b, 3) + 0, round(lo, 3) + 0, round(hi, 3) + 0)
fmt_sci  <- function(p) formatC(p, format = "E", digits = 1)
fmt_n    <- function(n) format(n, big.mark = ",")
fmt_pct  <- function(x) sprintf("%.0f%%", 100 * x)

n_obs_disc <- max(s5_num$N_disc, na.rm = TRUE)
n_obs_repl <- max(s5_num$N_repl, na.rm = TRUE)

hdr_disc <- paste0("Discovery set: Arizona + South Dakota (", fmt_n(n_obs_disc), " observations)")
hdr_repl <- paste0("Replication set: Oklahoma (", fmt_n(n_obs_repl), " observations)")
hdr_meta <- "Meta-analysis of both sets (fixed effect)"
hdr_het  <- "Heterogeneity between sets"

abbr_all <- c(AC = "acylcarnitine", CE = "cholesteryl ester", CER = "ceramide", DAG = "diacylglycerol",
              FA = "fatty acid", GlcCer = "glucosylceramide", LPC = "lysophosphatidylcholine",
              LPE = "lysophosphatidylethanolamine", PC = "phosphatidylcholine",
              PE = "phosphatidylethanolamine", PG = "phosphatidylglycerol", PI = "phosphatidylinositol",
              PS = "phosphatidylserine", SM = "sphingomyelin", TAG = "triacylglycerol")
cls_in   <- levels(droplevels(s5_num$Class_short))
abbr_txt <- paste0(names(abbr_all)[names(abbr_all) %in% cls_in], ", ",
                   abbr_all[names(abbr_all) %in% cls_in], collapse = "; ")

note_split <- paste0("Sensitivity analysis splitting the sample by study center. Lipid levels, UPF intake and ",
                     "continuous covariates were standardized within each set, and no family had members ",
                     "in more than one center.")
note_beta  <- paste0("\u03B2 and 95% confidence intervals (CI) are the difference in lipid level (SD units) ",
                     "per 1-SD higher UPF intake (% of total energy), from linear mixed models pooling the ",
                     "baseline (2001-2003) and follow-up (2006-2009) examinations, with random intercepts ",
                     "for family and participant.")
note_adj   <- paste0("Models were adjusted for age, sex, study center (discovery set only), education, ",
                     "body mass index, smoking, physical activity, alcohol intake, lipid-lowering medication ",
                     "use, total energy intake and examination.")
note_sel   <- paste0("The table lists the ", nrow(s5_num), " known lipids associated with UPF intake at ",
                     "P < 0.05 in the discovery set (out of ", fmt_n(n_known_tested), " known lipids tested; about ",
                     round(0.05 * n_known_tested), " expected by chance). ",
                     "Replicated: same direction of association and P < 0.05 in the replication set (",
                     sum(s5_num$replicated), " of ", nrow(s5_num), "). ",
                     "In the 42: one of the 42 known lipids associated with UPF intake at FDR q < 0.05 in the ",
                     "main analysis of all centers (Supplementary Table S4; ", sum(s5_num$in42), " of ",
                     nrow(s5_num), ").")
note_meta  <- paste0("Meta-analysis: fixed-effect inverse-variance weighted average of the discovery and ",
                     "replication estimates. Because the two sets together form the full sample, the meta-analyzed ",
                     "estimates are expected to be close to those of the main analysis (Supplementary Table S4).")
note_het   <- paste0("Heterogeneity: P-het is from Cochran's Q test (1 degree of freedom), which tests whether the ",
                     "discovery and replication estimates differ; I\u00B2 = (Q \u2212 1)/Q, set to 0% when Q \u2264 1, ",
                     "is the share of the difference between the two estimates beyond chance. ",
                     "Median I\u00B2 = ", fmt_pct(med_I2), "; ", n_het05, " of ", n_hits,
                     " lipids had P-het < 0.05 (about ", round(0.05 * n_hits, 1), " expected by chance).")
note_order <- paste0("Lipids with a positive discovery \u03B2 are listed first, then those with an inverse ",
                     "discovery \u03B2; within each group, lipids are ordered by class and then by name.")
note_abbr  <- paste0("Abbreviations: ", abbr_txt, "; UPF, ultra-processed food.")
notes_S5   <- c(note_split, note_beta, note_adj, note_sel, note_meta, note_het, note_order, note_abbr)

s5_fmt <- s5_num %>%
  transmute(Direction,
            class       = as.character(Class_short),
            lipid       = lipids_name,
            beta_disc   = fmt_beta(beta_disc, lower_CI_disc, upper_CI_disc),
            p_disc      = fmt_sci(p_disc),
            beta_repl   = fmt_beta(beta_repl, lower_CI_repl, upper_CI_repl),
            p_repl      = fmt_sci(p_repl),
            beta_meta   = fmt_beta(beta_meta, lower_CI_meta, upper_CI_meta),
            p_meta      = fmt_sci(p_meta),
            I2_het      = fmt_pct(I2_het),
            p_het       = fmt_sci(p_het),
            replicated  = ifelse(replicated, "Yes", "No"),
            in42        = ifelse(in42, "Yes", "No"))
col_names <- c("Lipid class", "Lipid", "\u03B2 (95% CI)", "P", "\u03B2 (95% CI)", "P",
               "\u03B2 (95% CI)", "P", "I\u00B2", "P-het",
               "Replicated", "In the 42 UPF-related lipids")
cols_disc <- 3:4; cols_repl <- 5:6; cols_meta <- 7:8; cols_het <- 9:10
cols_single <- c(1, 2, 11, 12)

n_pos   <- sum(s5_fmt$Direction == "Positive")
n_inv   <- sum(s5_fmt$Direction == "Inverse")
tab_pos <- s5_fmt %>% filter(Direction == "Positive") %>% select(-Direction)
tab_inv <- s5_fmt %>% filter(Direction == "Inverse")  %>% select(-Direction)
sub_pos <- paste0("Positively associated with UPF intake in the discovery set (n = ", n_pos, ")")
sub_inv <- paste0("Inversely associated with UPF intake in the discovery set (n = ", n_inv, ")")

wb <- createWorkbook()
addWorksheet(wb, "Table S5")
n_col <- length(col_names)

row_grp       <- 1
row_hdr       <- 2
row_sub_pos   <- 3
row_pos_start <- 4
row_pos_end   <- row_pos_start + n_pos - 1
row_sub_inv   <- row_pos_end + 1
row_inv_start <- row_sub_inv + 1
row_inv_end   <- row_inv_start + n_inv - 1
note_start    <- row_inv_end + 2

grp_list <- list(list(cols = cols_disc, txt = hdr_disc), list(cols = cols_repl, txt = hdr_repl),
                 list(cols = cols_meta, txt = hdr_meta), list(cols = cols_het,  txt = hdr_het))
for (g in grp_list) {
  writeData(wb, "Table S5", g$txt, startRow = row_grp, startCol = min(g$cols))
  mergeCells(wb, "Table S5", cols = g$cols, rows = row_grp)
}
writeData(wb, "Table S5", t(col_names[3:10]), startRow = row_hdr, startCol = 3, colNames = FALSE)
for (cc in cols_single) {
  writeData(wb, "Table S5", col_names[cc], startRow = row_grp, startCol = cc)
  mergeCells(wb, "Table S5", cols = cc, rows = row_grp:row_hdr)
}

writeData(wb, "Table S5", sub_pos, startRow = row_sub_pos, startCol = 1)
writeData(wb, "Table S5", tab_pos, startRow = row_pos_start, startCol = 1, colNames = FALSE)
writeData(wb, "Table S5", sub_inv, startRow = row_sub_inv, startCol = 1)
writeData(wb, "Table S5", tab_inv, startRow = row_inv_start, startCol = 1, colNames = FALSE)
writeData(wb, "Table S5", notes_S5, startRow = note_start, startCol = 1)

header_style <- createStyle(textDecoration = "bold", halign = "center", valign = "center", wrapText = TRUE)
sub_style    <- createStyle(textDecoration = c("bold", "italic"), halign = "left", valign = "center")
left_style   <- createStyle(halign = "left",   valign = "center")
num_style    <- createStyle(halign = "center", valign = "center")
note_style   <- createStyle(wrapText = TRUE, valign = "top", fontSize = 10)

addStyle(wb, "Table S5", header_style, rows = row_grp:row_hdr, cols = 1:n_col, gridExpand = TRUE)
addStyle(wb, "Table S5", createStyle(border = "top"), rows = row_grp, cols = 1:n_col,
         gridExpand = TRUE, stack = TRUE)
addStyle(wb, "Table S5", createStyle(border = "bottom", borderStyle = "hair"), rows = row_grp, cols = 3:10,
         gridExpand = TRUE, stack = TRUE)
addStyle(wb, "Table S5", createStyle(border = "bottom"), rows = row_hdr, cols = 1:n_col,
         gridExpand = TRUE, stack = TRUE)
for (r in c(row_sub_pos, row_sub_inv)) {
  mergeCells(wb, "Table S5", cols = 1:n_col, rows = r)
  addStyle(wb, "Table S5", sub_style, rows = r, cols = 1:n_col, gridExpand = TRUE)
}
body_rows <- c(row_pos_start:row_pos_end, row_inv_start:row_inv_end)
addStyle(wb, "Table S5", left_style, rows = body_rows, cols = 1:2,     gridExpand = TRUE)
addStyle(wb, "Table S5", num_style,  rows = body_rows, cols = 3:n_col, gridExpand = TRUE)
addStyle(wb, "Table S5", createStyle(border = "bottom"), rows = row_inv_end, cols = 1:n_col,
         gridExpand = TRUE, stack = TRUE)

for (r in note_start:(note_start + length(notes_S5) - 1)) {
  mergeCells(wb, "Table S5", cols = 1:n_col, rows = r)
  addStyle(wb, "Table S5", note_style, rows = r, cols = 1:n_col, gridExpand = TRUE)
}
setRowHeights(wb, "Table S5", rows = row_grp, heights = 30)
setRowHeights(wb, "Table S5", rows = note_start:(note_start + length(notes_S5) - 1), heights = 45)

setColWidths(wb, "Table S5", cols = 1:n_col, widths = c(11, 18, 26, 10, 26, 10, 26, 10, 8, 10, 11, 14))
freezePane(wb, "Table S5", firstActiveRow = row_sub_pos)

setwd(tabsFigs)
saveWorkbook(wb, "tableS5_sensCenter_OK.xlsx", overwrite = TRUE)
setwd(result)
write.csv(s5_num, "tableS5_sensCenter_OK_numbers.csv", row.names = FALSE)
