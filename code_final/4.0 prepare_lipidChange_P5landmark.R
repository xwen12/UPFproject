# Build the follow-up (Phase 5) landmark dataset for 5-year lipid change analyses (4.1-4.3, 5.6).
# Input: analysis datasets (1.3), 2.1 results. Output: lipid change dataset (data/).

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
raw_dir  <- "path/to/SHS_raw_data"
library(dplyr)
library(lubridate)
rawData <- raw_dir
readyData <- file.path(proj_dir, "data")
result <-file.path(proj_dir, "result")

setwd(readyData)
load("lipidomics_p4.Rda")
load("lipidomics_p5.Rda")
load("dat_AimOne_p4.Rda")
load("dat_AimOne_p5.Rda")
setwd(result)
key_lipids_risk <- read.table('UPF_lipids_risk.csv',sep=',',header=T)
key_lipids_prot <- read.table('UPF_lipids_prot.csv',sep=',',header=T)
key_lipids <- c(key_lipids_risk$lipids_name, key_lipids_prot$lipids_name)
setwd(rawData)
cvd_raw <- read.csv("cvd data_2020.csv", header = TRUE, stringsAsFactors = FALSE)
t2d_p7  <- read.csv("p7data/p7_analytic_dataset.csv", header = TRUE, stringsAsFactors = FALSE)

lip4 <- as.data.frame(lipidomics_p4)[, c("ID", key_lipids)]
lip5 <- as.data.frame(lipidomics_p5)[, c("ID", key_lipids)]
lip4[key_lipids] <- lapply(lip4[key_lipids], as.numeric)
lip5[key_lipids] <- lapply(lip5[key_lipids], as.numeric)

lip45 <- inner_join(lip4, lip5, by = "ID", suffix = c("_p4", "_p5"))

lip_change <- data.frame(ID = lip45$ID)
for(lp in key_lipids){
  lip_change[[paste0(lp, "_bsl")]] <- lip45[[paste0(lp, "_p4")]]
  lip_change[[paste0(lp, "_chg")]] <- lip45[[paste0(lp, "_p5")]] - lip45[[paste0(lp, "_p4")]]
}

batch_check <- data.frame(
  lipids_name = key_lipids,
  mean_P4     = sapply(key_lipids, function(lp) mean(lip45[[paste0(lp, "_p4")]], na.rm = TRUE)),
  mean_change = sapply(key_lipids, function(lp) mean(lip_change[[paste0(lp, "_chg")]], na.rm = TRUE)),
  sd_change   = sapply(key_lipids, function(lp) sd(lip_change[[paste0(lp, "_chg")]],   na.rm = TRUE)),
  r_p4_p5     = sapply(key_lipids, function(lp) cor(lip45[[paste0(lp, "_p4")]], lip45[[paste0(lp, "_p5")]],
                                                    use = "complete.obs", method = "spearman")))
batch_check$pct_change <- 100 * batch_check$mean_change / batch_check$mean_P4
setwd(result)
write.csv(batch_check, "lipidChange_batchCheck.csv", row.names = FALSE)

cov5 <- dat_AimOne_p5 %>%
  select(IDNO, S5AGE, S5EDU, S5BMI, S5SMOKE, S5ETOH, ped_avg3day,
         LIPIDLOW_p5 = LIPIDLOW, S5ADADM2020, S5USHTN, S5DT_KCAL, S5G0)
cov4 <- dat_AimOne_p4 %>%
  select(IDNO, FID, SEX, CENTER, S4BMI, S4DT_KCAL, S4G0, S4ADADM2020)
dat_1 <- cov4 %>%
  inner_join(cov5, by = "IDNO") %>%
  inner_join(lip_change, by = c("IDNO" = "ID")) %>%
  mutate(BMI_diff     = S5BMI     - S4BMI,
         DT_KCAL_diff = S5DT_KCAL - S4DT_KCAL,
         G0_diff      = S5G0      - S4G0)

parse_shs_date <- function(x){
  x <- as.character(x)
  x[x == ""] <- NA
  as.Date(lubridate::parse_date_time(x, orders = c("mdY", "Ymd", "dbY"), quiet = TRUE))
}
cvd5 <- cvd_raw %>%
  select(IDNO, S5EXAMDT, DEADCODE, DOD,
         S5CVDFREE, ANYCVDS5, ANYCVDDTS5,
         S5CHDFREE, ANYCHDS5, ANYCHDDTS5,
         S5MIFREE,  ANYMIS5,  ANYMIDTS5,
         S5DEFSTK,  ANYDSTKS5, ANYDSTKDTS5) %>%
  mutate(across(c(S5EXAMDT, DOD, ANYCVDDTS5, ANYCHDDTS5, ANYMIDTS5, ANYDSTKDTS5), parse_shs_date),
         DEADCODE = as.character(DEADCODE))
for(v in c("S5EXAMDT","DOD","ANYCVDDTS5","ANYCHDDTS5","ANYMIDTS5","ANYDSTKDTS5")){
  raw_v <- as.character(cvd_raw[[v]])
}

closing_date <- as.Date("2020-12-31")
dat_2 <- dat_1 %>%
  left_join(cvd5, by = "IDNO") %>%
  mutate(end_fu = pmin(DOD, closing_date, na.rm = TRUE))

make_outcome <- function(event, event_date, dod, end_fu, s5date){
  date_use <- case_when(event == 1 & !is.na(event_date) ~ event_date,
                        event == 1 &  is.na(event_date) ~ dod,
                        TRUE                             ~ end_fu)
  as.numeric(date_use - s5date)
}
dat_2 <- dat_2 %>%
  mutate(
    chd_death_after5 = !is.na(DEADCODE) & DEADCODE %in% c("2","02") & !is.na(DOD) & DOD > S5EXAMDT,
    atrisk_chd5 = !is.na(S5CHDFREE) & S5CHDFREE == 1 & !is.na(S5MIFREE) & S5MIFREE == 1,
    chd5        = as.integer((ANYCHDS5 %in% 1) | (ANYMIS5 %in% 1) | chd_death_after5),
    chd5_date   = pmin(ANYCHDDTS5, ANYMIDTS5, na.rm = TRUE),
    TIME_chd5   = make_outcome(chd5, chd5_date, DOD, end_fu, S5EXAMDT),
    atrisk_cvd5 = !is.na(S5CVDFREE) & S5CVDFREE == 1 & atrisk_chd5,
    cvd5        = as.integer((ANYCVDS5 %in% 1) | chd5 == 1),
    cvd5_date   = pmin(ANYCVDDTS5, chd5_date, na.rm = TRUE),
    TIME_cvd5   = make_outcome(cvd5, cvd5_date, DOD, end_fu, S5EXAMDT),
    atrisk_mi5  = !is.na(S5MIFREE) & S5MIFREE == 1,
    mi5         = as.integer(ANYMIS5 %in% 1),
    TIME_mi5    = make_outcome(mi5, ANYMIDTS5, DOD, end_fu, S5EXAMDT),
    atrisk_stroke5 = !is.na(S5DEFSTK) & S5DEFSTK == 0,
    stroke5        = as.integer(ANYDSTKS5 %in% 1),
    TIME_stroke5   = make_outcome(stroke5, ANYDSTKDTS5, DOD, end_fu, S5EXAMDT)
  )

load(file.path(readyData, "outcome_mi_stroke_p4.Rda"))
dat_2 <- dat_2 %>%
  left_join(dat_AimOne_p4 %>% select(IDNO, cvd, chd, TIME_cvd, TIME_chd), by = "IDNO") %>%
  left_join(outcome_mi_stroke_p4, by = "IDNO") %>%
  left_join(cvd_raw %>% select(IDNO, S4EXAMDT) %>% mutate(S4EXAMDT = parse_shs_date(S4EXAMDT)), by = "IDNO") %>%
  mutate(day_p5 = as.numeric(S5EXAMDT - S4EXAMDT),
         atrisk_cvd5    = atrisk_cvd5    & !(cvd %in% TRUE & as.numeric(TIME_cvd) <= day_p5),
         atrisk_chd5    = atrisk_chd5    & !(chd %in% TRUE & as.numeric(TIME_chd) <= day_p5),
         atrisk_mi5     = atrisk_mi5     & !(mi  %in% 1    & TIME_mi     <= day_p5),
         atrisk_stroke5 = atrisk_stroke5 & !(definite_stroke %in% 1 & TIME_stroke <= day_p5))

for(o in c("cvd5","chd5","mi5","stroke5")){
  ar <- dat_2[[paste0("atrisk_", o)]]; tt <- dat_2[[paste0("TIME_", o)]]
}

t2d_p7 <- t2d_p7 %>%
  select(IDNO, S7ADADM2025) %>%
  filter(!is.na(S7ADADM2025), S7ADADM2025 != "")
dat_3 <- dat_2 %>%
  left_join(t2d_p7, by = "IDNO") %>%
  mutate(atrisk_dm7 = !is.na(S5ADADM2020) & S5ADADM2020 != "DM" &
                      !is.na(S4ADADM2020) & S4ADADM2020 != "DM" &
                      !is.na(S7ADADM2025),
         incidentDM_p5p7 = ifelse(S7ADADM2025 == "DM", 1L, 0L))

dat_lipidChange <- dat_3 %>%
  select(IDNO, FID, SEX, CENTER,
         S5AGE, S5EDU, S5BMI, S5SMOKE, S5ETOH, ped_avg3day, LIPIDLOW_p5,
         S5ADADM2020, S5USHTN, S5DT_KCAL, S5G0, S4G0,
         BMI_diff, DT_KCAL_diff, G0_diff,
         atrisk_cvd5, cvd5, TIME_cvd5, atrisk_chd5, chd5, TIME_chd5,
         atrisk_mi5, mi5, TIME_mi5, atrisk_stroke5, stroke5, TIME_stroke5,
         atrisk_dm7, incidentDM_p5p7,
         ends_with("_bsl"), ends_with("_chg"))
setwd(readyData)
save(dat_lipidChange, file = "dat_lipidChange_P5landmark.Rda")
