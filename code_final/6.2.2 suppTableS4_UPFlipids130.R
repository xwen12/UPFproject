# Supplemental Table S4: the 130 lipids associated with UPF intake at FDR < 0.05 (42 named + 88 unknown).
# Input: 2.1 results, lipid class annotation. Output: Table S4 (result/tabsFigs).

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
P4P5all <- read.table("P4P5all.csv", sep = ",", header = TRUE, check.names = FALSE)

upf_lipids42 <- read.table("UPF_lipids_all.csv", sep = ",", header = TRUE, check.names = FALSE)

s4_raw <- P4P5all %>%
  filter(qvals < 0.05) %>%
  mutate(known = as.logical(known))

lipid_info <- read.csv(lipid_info_path, header = TRUE, stringsAsFactors = FALSE,
                       check.names = FALSE)

lipid_info <- lipid_info %>%
  select(lipids_name, Class_short) %>%
  distinct(lipids_name, .keep_all = TRUE)

s4_dt <- s4_raw %>%
  left_join(lipid_info, by = "lipids_name")

s4_dt <- s4_dt %>%
  mutate(Class_short = ifelse(is.na(Class_short) | Class_short == "", "N/A", Class_short))
s4_dt[s4_dt$Class_short == "N/A",                       "Class_short"] <- "Unknown"
s4_dt[s4_dt$Class_short %in% c("FAHFA", "scFA", "USFA"), "Class_short"] <- "FA"
s4_dt[s4_dt$Class_short %in% c("Glc"),                  "Class_short"] <- "GlcCer"

cls_order <- c("AC", "CE", "CER", "DAG", "FA", "GlcCer", "PC", "PE", "PG", "PI",
               "PS", "SM", "TAG", "Unknown")
cls_extra <- sort(setdiff(unique(s4_dt$Class_short), cls_order))
cls_order <- c(setdiff(cls_order, "Unknown"), cls_extra, "Unknown")

s4_num <- s4_dt %>%
  mutate(Direction   = factor(ifelse(beta > 0, "Positive", "Inverse"),
                              levels = c("Positive", "Inverse")),
         Class_short = factor(Class_short, levels = cls_order),
         is_unk_id   = str_detect(lipids_name, "^Unknown_[0-9]+$"),
         name_rank   = order(str_order(lipids_name, numeric = TRUE))) %>%
  arrange(Direction, Class_short, is_unk_id, name_rank) %>%
  select(Direction, Class_short, lipids_name, known, beta, se, lower_CI, upper_CI, pvals, qvals)

fmt_beta <- function(b, lo, hi) sprintf("%.3f (%.3f, %.3f)", b, lo, hi)
fmt_sci  <- function(p) formatC(p, format = "E", digits = 1)

note_beta  <- paste0("\u03B2 and 95% confidence intervals (CI) are the difference in lipid level (SD units) ",
                     "per 1-SD higher UPF intake (% of total energy), from linear mixed models pooling the ",
                     "baseline (2001-2003) and follow-up (2006-2009) examinations, with random intercepts ",
                     "for family and participant.")
note_adj   <- paste0("Models were adjusted for age, sex, study center, education, body mass index, smoking, ",
                     "physical activity, alcohol intake, lipid-lowering medication use, total energy intake ",
                     "and examination.")
note_fdr   <- paste0("FDR q: Benjamini-Hochberg false discovery rate across all 1,542 lipids tested. ",
                     "The table lists all ", nrow(s4_num), " lipids with FDR q < 0.05 (",
                     sum(s4_num$known), " known and ", sum(!s4_num$known), " unknown).")
note_order <- paste0("Positive markers (\u03B2 > 0) are listed first, then inverse markers (\u03B2 < 0); ",
                     "within each group, lipids are ordered by class and then by name. ",
                     "Unknown: lipid features without a confirmed class annotation.")
note_abbr  <- paste0("Abbreviations: AC, acylcarnitine; CE, cholesteryl ester; CER, ceramide; ",
                     "DAG, diacylglycerol; FA, fatty acid; GlcCer, glucosylceramide; ",
                     "PC, phosphatidylcholine; PE, phosphatidylethanolamine; PG, phosphatidylglycerol; ",
                     "PI, phosphatidylinositol; PS, phosphatidylserine; SM, sphingomyelin; ",
                     "TAG, triacylglycerol; UPF, ultra-processed food.")
notes_S4   <- c(note_beta, note_adj, note_fdr, note_order, note_abbr)

s4_fmt <- s4_num %>%
  transmute(Direction,
            `Lipid class`         = as.character(Class_short),
            `Lipid`               = lipids_name,
            beta_ci               = fmt_beta(beta, lower_CI, upper_CI),
            `P`                   = fmt_sci(pvals),
            `FDR q`               = fmt_sci(qvals))
names(s4_fmt)[names(s4_fmt) == "beta_ci"] <- "\u03B2 (95% CI)"

n_pos   <- sum(s4_fmt$Direction == "Positive")
n_inv   <- sum(s4_fmt$Direction == "Inverse")
tab_pos <- s4_fmt %>% filter(Direction == "Positive") %>% select(-Direction)
tab_inv <- s4_fmt %>% filter(Direction == "Inverse")  %>% select(-Direction)
sub_pos <- paste0("Positively associated with UPF intake (n = ", n_pos, ")")
sub_inv <- paste0("Inversely associated with UPF intake (n = ", n_inv, ")")

wb <- createWorkbook()
addWorksheet(wb, "Table S4")
n_col <- ncol(tab_pos)

row_sub_pos   <- 2
row_pos_start <- 3
row_pos_end   <- row_pos_start + n_pos - 1
row_sub_inv   <- row_pos_end + 1
row_inv_start <- row_sub_inv + 1
row_inv_end   <- row_inv_start + n_inv - 1
note_start    <- row_inv_end + 2

writeData(wb, "Table S4", t(names(tab_pos)), startRow = 1, startCol = 1, colNames = FALSE)
writeData(wb, "Table S4", sub_pos, startRow = row_sub_pos, startCol = 1)
writeData(wb, "Table S4", tab_pos, startRow = row_pos_start, startCol = 1, colNames = FALSE)
writeData(wb, "Table S4", sub_inv, startRow = row_sub_inv, startCol = 1)
writeData(wb, "Table S4", tab_inv, startRow = row_inv_start, startCol = 1, colNames = FALSE)
writeData(wb, "Table S4", notes_S4, startRow = note_start, startCol = 1)

header_style <- createStyle(textDecoration = "bold", halign = "center", valign = "center",
                            border = "TopBottom", wrapText = TRUE)
sub_style    <- createStyle(textDecoration = c("bold", "italic"), halign = "left", valign = "center")
left_style   <- createStyle(halign = "left",   valign = "center")
num_style    <- createStyle(halign = "center", valign = "center")
note_style   <- createStyle(wrapText = TRUE, valign = "top", fontSize = 10)

addStyle(wb, "Table S4", header_style, rows = 1, cols = 1:n_col, gridExpand = TRUE)
for (r in c(row_sub_pos, row_sub_inv)) {
  mergeCells(wb, "Table S4", cols = 1:n_col, rows = r)
  addStyle(wb, "Table S4", sub_style, rows = r, cols = 1:n_col, gridExpand = TRUE)
}
body_rows <- c(row_pos_start:row_pos_end, row_inv_start:row_inv_end)
addStyle(wb, "Table S4", left_style, rows = body_rows, cols = 1:2,     gridExpand = TRUE)
addStyle(wb, "Table S4", num_style,  rows = body_rows, cols = 3:n_col, gridExpand = TRUE)
addStyle(wb, "Table S4", createStyle(border = "bottom"), rows = row_inv_end, cols = 1:n_col,
         gridExpand = TRUE, stack = TRUE)

for (r in note_start:(note_start + length(notes_S4) - 1)) {
  mergeCells(wb, "Table S4", cols = 1:n_col, rows = r)
  addStyle(wb, "Table S4", note_style, rows = r, cols = 1:n_col, gridExpand = TRUE)
}
setRowHeights(wb, "Table S4", rows = note_start:(note_start + length(notes_S4) - 1), heights = 45)

setColWidths(wb, "Table S4", cols = 1:n_col, widths = c(12, 24, 26, 11, 11))
freezePane(wb, "Table S4", firstRow = TRUE)

setwd(tabsFigs)
saveWorkbook(wb, "tableS4_UPFlipids130.xlsx", overwrite = TRUE)
setwd(result)
write.csv(s4_num, "tableS4_UPFlipids130_numbers.csv", row.names = FALSE)
