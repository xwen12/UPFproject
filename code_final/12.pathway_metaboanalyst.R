# Prepare the 42 UPF-related lipids for KEGG pathway analysis (MetaboAnalyst 6.0) and LIPEA.
# Input: result/UPF_lipids_all.csv (2.1) + HMDB ID map. Output: lipid lists for upload.

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
raw_dir  <- "path/to/SHS_raw_data"
library(dplyr)
library(stringr)

rawData       <- raw_dir
result        <- file.path(proj_dir, "result")
enrich_output <- file.path(proj_dir, "result", "enrich")
hmdb_map_path <- file.path(raw_dir, "lipidomics data", "lipids_IDs_XiaoxiaoEdit_filled.csv")

kegg_col <- NULL

if (!dir.exists(enrich_output)) dir.create(enrich_output, recursive = TRUE)

setwd(result)

P4P5all <- read.table('P4P5all.csv', sep = ',', header = TRUE)

setA <- read.table('UPF_lipids_all.csv', sep = ',', header = TRUE)

setB <- P4P5all %>%
  filter(known, pvals < 0.05)

bg <- P4P5all %>%
  filter(known)

hmdb_map <- read.csv(hmdb_map_path, header = TRUE, stringsAsFactors = FALSE,
                     na.strings = c("", "NA", "N/A"))

if (is.null(kegg_col)) {
  kegg_col <- grep("kegg", names(hmdb_map), ignore.case = TRUE, value = TRUE)
}
hmdb_map <- hmdb_map %>%
  rename(KEGG.ID = all_of(kegg_col)) %>%
  mutate(KEGG.ID = as.character(KEGG.ID))

hmdb_map_clean <- hmdb_map %>%
  distinct(lipids_me, .keep_all = TRUE) %>%
  select(lipids_me, HMDB.ID, KEGG.ID, Pubchem.ID) %>%
  mutate(HMDB_multi = str_detect(HMDB.ID, "[;,|]") & !is.na(HMDB.ID),
         HMDB.ID    = str_trim(str_split_fixed(HMDB.ID, "[;,|]", 2)[, 1]),
         HMDB.ID    = na_if(HMDB.ID, "")) %>%
  mutate(HMDB.ID = ifelse(
    str_detect(HMDB.ID, "^HMDB[0-9]+$"),
    sprintf("HMDB%07d", as.integer(str_remove(HMDB.ID, "^HMDB"))),
    HMDB.ID))

hmdb_map_clean <- hmdb_map_clean %>%
  mutate(KEGG_multi = str_detect(KEGG.ID, "[;,|]") & !is.na(KEGG.ID),
         KEGG.ID    = str_trim(str_split_fixed(KEGG.ID, "[;,|]", 2)[, 1]),
         KEGG.ID    = str_remove(KEGG.ID, "^(cpd:|CPD:)"),
         KEGG.ID    = na_if(KEGG.ID, ""),
         KEGG_bad   = !is.na(KEGG.ID) & !str_detect(KEGG.ID, "^C[0-9]{5}$"))

hmdb_map_clean <- hmdb_map_clean %>%
  mutate(KEGG.ID = ifelse(KEGG_bad, NA, KEGG.ID))

strip_isomer <- function(nm) {
  gsub("\\)\\s*[A-Z]$", ")", trimws(nm), perl = TRUE)
}

audit <- bg %>%
  select(lipids_name, beta, pvals, qvals) %>%
  mutate(lipid_clean = strip_isomer(lipids_name),
         in_setA     = lipids_name %in% setA$lipids_name,
         in_setB     = lipids_name %in% setB$lipids_name) %>%
  left_join(hmdb_map_clean, by = c("lipids_name" = "lipids_me")) %>%
  mutate(in_map   = lipids_name %in% hmdb_map_clean$lipids_me,
         has_HMDB = !is.na(HMDB.ID),
         has_KEGG = !is.na(KEGG.ID))

setwd(enrich_output)
write.table(audit, 'lipid_HMDB_mapping_audit.csv', sep = ',', row.names = F, col.names = T)

hmdb_A  <- audit %>% filter(in_setA, has_HMDB) %>% pull(HMDB.ID) %>% unique()
hmdb_B  <- audit %>% filter(in_setB, has_HMDB) %>% pull(HMDB.ID) %>% unique()
hmdb_bg <- audit %>% filter(has_HMDB)          %>% pull(HMDB.ID) %>% unique()

save_txt <- function(vec, filename) {
  writeLines(vec, filename)
}
save_txt(hmdb_A,  'MA_target_HMDB_q005.txt')
save_txt(hmdb_B,  'MA_target_HMDB_p005.txt')
save_txt(hmdb_bg, 'MA_background_HMDB.txt')

kegg_A  <- audit %>% filter(in_setA, has_KEGG) %>% pull(KEGG.ID) %>% unique()
kegg_B  <- audit %>% filter(in_setB, has_KEGG) %>% pull(KEGG.ID) %>% unique()
kegg_bg <- audit %>% filter(has_KEGG)          %>% pull(KEGG.ID) %>% unique()

save_txt(kegg_A,  'MA_target_KEGG_q005.txt')
save_txt(kegg_B,  'MA_target_KEGG_p005.txt')
save_txt(kegg_bg, 'MA_background_KEGG.txt')

save_name <- function(vec, filename) {
  write.csv(data.frame(lipids_name = vec), file = filename,
            row.names = FALSE, quote = TRUE)
}
save_name(audit$lipid_clean[audit$in_setA], 'LION_target_name_q005.csv')
save_name(audit$lipid_clean[audit$in_setB], 'LION_target_name_p005.csv')
save_name(audit$lipid_clean,                'LION_background_name.csv')
