# Count named lipids per class and UPF-related lipids per class for the Figure 1 donut chart.
# Input: result/P4P5all.csv (2.1), lipid class annotation. Output: result/fig1_lipidClassCounts.csv.

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
raw_dir  <- "path/to/SHS_raw_data"

library(dplyr)
result <- file.path(proj_dir, "result")
lipid_info_path <- file.path(raw_dir, "lipidomics data", "lipid_4_6_compare.csv")

setwd(result)
P4P5all    <- read.csv("P4P5all.csv", check.names = FALSE)
lipid_info <- read.csv(lipid_info_path, stringsAsFactors = FALSE, check.names = FALSE)

dt0 <- P4P5all %>%
  filter(known) %>%
  select(lipids_name, qvals, beta) %>%
  left_join(lipid_info[, c("lipids_name", "Class_short")], by = "lipids_name") %>%
  mutate(Class_short = case_when(
    is.na(Class_short) | Class_short == "N/A" ~ "N/A",
    Class_short %in% c("FAHFA", "scFA", "USFA") ~ "FA",
    Class_short == "Glc"                        ~ "GlcCer",
    TRUE                                        ~ Class_short))

cls_n <- dt0 %>%
  group_by(Class_short) %>%
  summarise(n_lipids = n(),
            n_upf42  = sum(qvals < 0.05), .groups = "drop") %>%
  arrange(desc(n_lipids))

write.csv(cls_n, "fig1_lipidClassCounts.csv", row.names = FALSE)
