# Supplemental Table S1: FFQ food items in each UPF subclass.
# Input: NOVA-labelled FFQ food lists (data/NOVAlabels). Output: Table S1 (result/tabsFigs).

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
library(tidyverse)
library(readxl)
library(openxlsx)

readyData <- file.path(proj_dir, "data")
tabsFigs  <- file.path(proj_dir, "result", "tabsFigs")

if (!dir.exists(tabsFigs)) dir.create(tabsFigs, recursive = TRUE)

setwd(readyData)
labels_p4 <- read_excel("NOVAlabels/p4_food_list_final.xlsx", sheet = "P4")
labels_p5 <- read_excel("NOVAlabels/p5_food_list_final.xlsx", sheet = "Sheet1")

names(labels_p4) <- trimws(names(labels_p4))
names(labels_p5) <- trimws(names(labels_p5))

labels_all <- bind_rows(
  labels_p4 %>% select(Variable, FoodLabel, NOVA_class, NOVA_Subclass) %>%
    mutate(Visit = "P4"),
  labels_p5 %>% select(Variable, FoodLabel, NOVA_class, NOVA_Subclass) %>%
    mutate(Visit = "P5")
)

subclass_names <- c(
  "1" = "Ultra-processed breads & breakfast foods",
  "2" = "Sauces, spreads & condiments",
  "3" = "Packaged sweet snacks & desserts",
  "4" = "Packaged savory snacks",
  "5" = "Sugar/artificially sweetened beverages",
  "6" = "Meat/poultry/seafood-based ready-to-eat products",
  "7" = "Ready-to-eat/heat mixed dishes",
  "8" = "Dairy-based desserts"
)

subclass_order <- c(
  "Ready-to-eat/heat mixed dishes",
  "Packaged sweet snacks & desserts",
  "Sugar/artificially sweetened beverages",
  "Sauces, spreads & condiments",
  "Meat/poultry/seafood-based ready-to-eat products",
  "Ultra-processed breads & breakfast foods",
  "Packaged savory snacks",
  "Dairy-based desserts"
)

liquor_labels <- c("Liquor or mixed drinks (gm)", "Liquor")

note_ffq     <- paste0("Food items are shown as labeled in the food frequency questionnaire (FFQ) used at each visit. ",
                       "The FFQ differed between baseline (2001-2003) and follow-up (2006-2009), ",
                       "so the item lists are not the same.")
note_subcls  <- paste0("UPF was defined as NOVA group 4. UPF subclasses were adapted from Fang Z, et al. ",
                       "BMJ 2024;385:e078476.")
note_liquor  <- paste0("Liquor or mixed drinks (baseline) and liquor (follow-up) were classified as NOVA group 4 ",
                       "but were not included in total UPF or UPF subclass intake.")
note_abbr    <- "Abbreviations: FFQ, food frequency questionnaire; UPF, ultra-processed food."

upf_items <- labels_all %>%
  filter(NOVA_class == 4) %>%
  filter(!FoodLabel %in% liquor_labels) %>%
  mutate(
    FoodLabel = str_remove(FoodLabel, "\\s*\\(gm\\)\\s*$"),
    FoodLabel = str_squish(FoodLabel),
    Subclass  = subclass_names[as.character(NOVA_Subclass)]
  ) %>%
  arrange(Visit, Variable)

s1_long <- upf_items %>%
  group_by(Subclass, Visit) %>%
  summarise(items   = paste(FoodLabel, collapse = "; "),
            .groups = "drop")

tabS1 <- s1_long %>%
  pivot_wider(names_from = Visit, values_from = items) %>%
  mutate(Subclass = factor(Subclass, levels = subclass_order)) %>%
  arrange(Subclass) %>%
  transmute(
    `UPF subclass`                = as.character(Subclass),
    `Baseline FFQ food items`     = P4,
    `Follow-up FFQ food items`    = P5
  )

notes_S1 <- c(note_ffq, note_subcls, note_liquor, note_abbr)

wb <- createWorkbook()
addWorksheet(wb, "Table S1")

writeData(wb, "Table S1", tabS1, startRow = 1, startCol = 1)

note_start <- nrow(tabS1) + 3
writeData(wb, "Table S1", notes_S1, startRow = note_start, startCol = 1)

header_style <- createStyle(textDecoration = "bold", wrapText = TRUE, valign = "top",
                            border = "bottom")
body_style   <- createStyle(wrapText = TRUE, valign = "top")
note_style   <- createStyle(wrapText = TRUE, valign = "top", fontSize = 10)
addStyle(wb, "Table S1", header_style, rows = 1, cols = 1:3, gridExpand = TRUE)
addStyle(wb, "Table S1", body_style, rows = 2:(nrow(tabS1) + 1), cols = 1:3, gridExpand = TRUE)

for (r in note_start:(note_start + length(notes_S1) - 1)) {
  mergeCells(wb, "Table S1", cols = 1:3, rows = r)
  addStyle(wb, "Table S1", note_style, rows = r, cols = 1:3, gridExpand = TRUE)
}
setRowHeights(wb, "Table S1", rows = note_start:(note_start + length(notes_S1) - 1), heights = 30)

setColWidths(wb, "Table S1", cols = 1:3, widths = c(30, 80, 80))

setwd(tabsFigs)
saveWorkbook(wb, "tableS1_UPFsubclass_foodItems.xlsx", overwrite = TRUE)
