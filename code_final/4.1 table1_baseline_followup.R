# Table 1: participant characteristics at baseline and follow-up.
# Input: analysis datasets (1.3). Output: Table 1 (result/tabsFigs).

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
raw_dir  <- "path/to/SHS_raw_data"

library(dplyr)
library(tidyr)
library(tibble)

library(lmerTest)
library(geepack)
library(multgee)

library(flextable)
library(officer)

rawData   <- raw_dir
readyData <- file.path(proj_dir, "data")
result    <- file.path(proj_dir, "result")
tabsFigs  <- file.path(proj_dir, "result", "tabsFigs")

setwd(readyData)
load("lipidomics_p4.Rda")
load("lipidomics_p5.Rda")
load("dat_AimOne_p4.Rda")
load("dat_AimOne_p5.Rda")

setwd(rawData)
chd_raw <- read.csv("cvd data_2020.csv", header = TRUE, stringsAsFactors = FALSE)

ids_lip_p4 <- as.data.frame(lipidomics_p4)$ID
ids_lip_p5 <- as.data.frame(lipidomics_p5)$ID

p4 <- dat_AimOne_p4 %>% filter(IDNO %in% ids_lip_p4)
p5 <- dat_AimOne_p5 %>% filter(IDNO %in% ids_lip_p5)

cvd_subtypes <- chd_raw %>%
  select(IDNO, ANYDSTKS4, ANYMIS4) %>%
  mutate(ANYDSTKS4 = as.integer(ANYDSTKS4),
         ANYMIS4   = as.integer(ANYMIS4))

p4 <- p4 %>%
  left_join(cvd_subtypes, by = "IDNO") %>%
  mutate(
    definite_stroke = case_when(
      cvd == 1 & ANYDSTKS4 == 1 ~ 1L,
      cvd == 1 & ANYDSTKS4 != 1 ~ 0L,
      cvd == 0                  ~ 0L,
      TRUE                      ~ NA_integer_),
    mi = case_when(
      cvd == 1 & ANYMIS4 == 1 ~ 1L,
      cvd == 1 & ANYMIS4 != 1 ~ 0L,
      cvd == 0                ~ 0L,
      TRUE                    ~ NA_integer_)
  )

dat_b <- p4 %>%
  transmute(ID       = IDNO,
            FID      = as.character(FID),
            Phase    = "Baseline",
            AGE      = S4AGE,
            SEX, CENTER,
            EDU      = S4EDU,
            SMOKE    = S4SMOKE,
            ETOH     = S4ETOH,
            PA       = PDM3,
            DT_KCAL  = S4DT_KCAL,
            BMI      = S4BMI,
            SBP      = S4SBP,
            DBP      = S4DBP,
            USHTN    = S4USHTN,
            G0       = S4G0,
            ADADM    = S4ADADM2020,
            LIPIDLOW,
            UPF_pct  = nova4_kcal_P * 100)

dat_f <- p5 %>%
  transmute(ID       = IDNO,
            FID      = as.character(FID),
            Phase    = "Follow-up",
            AGE      = S5AGE,
            SEX, CENTER,
            EDU      = S5EDU,
            SMOKE    = S5SMOKE,
            ETOH     = S5ETOH,
            PA       = ped_avg3day,
            DT_KCAL  = S5DT_KCAL,
            BMI      = S5BMI,
            SBP      = S5SBP,
            DBP      = S5DBP,
            USHTN    = S5USHTN,
            G0       = S5G0,
            ADADM    = S5ADADM2020,
            LIPIDLOW,
            UPF_pct  = nova4_kcal_P * 100)

dat <- bind_rows(
  dat_b %>% mutate(across(where(is.factor), as.character)),
  dat_f %>% mutate(across(where(is.factor), as.character))
)

dat <- dat %>%
  mutate(
    ID       = as.factor(ID),
    FID      = as.factor(FID),
    Phase    = factor(Phase, levels = c("Baseline", "Follow-up")),
    SEX      = factor(SEX,      levels = c("M", "F"),
                      labels = c("Male", "Female")),
    CENTER   = factor(CENTER,   levels = c("AZ", "OK", "SD"),
                      labels = c("Arizona", "Oklahoma", "North and South Dakota")),
    SMOKE    = factor(SMOKE,    levels = c("N", "E", "Y"),
                      labels = c("Never", "Former", "Current")),
    ETOH     = factor(ETOH,     levels = c("N", "E", "Y"),
                      labels = c("Never", "Former", "Current")),
    USHTN    = factor(USHTN,    levels = c("N", "Y"),
                      labels = c("No", "Yes")),
    ADADM    = factor(ADADM,    levels = c("NFG", "IFG", "DM"),
                      labels = c("Normal fasting glucose",
                                 "Impaired fasting glucose",
                                 "Diabetes")),
    LIPIDLOW = factor(LIPIDLOW, levels = c("0", "1"),
                      labels = c("No", "Yes"))
  )

var_dict <- tribble(
  ~var,       ~label,                                  ~type,    ~level,   ~digits, ~test,
  "AGE",      "Age, years",                            "normal", NA,       1,       TRUE,
  "SEX",      "Female, n (%)",                         "binary", "Female", 1,       TRUE,
  "CENTER",   "Study center, n (%)",                   "multi",  NA,       1,       TRUE,
  "EDU",      "Education, years",                      "normal", NA,       1,       TRUE,
  "BMI",      "BMI, kg/m²",                       "normal", NA,       1,       TRUE,
  "SMOKE",    "Smoking status, n (%)",                 "multi",  NA,       1,       TRUE,
  "ETOH",     "Drinking status, n (%)",                "multi",  NA,       1,       TRUE,
  "PA",       "Physical activity, steps/day",          "normal", NA,       1,       TRUE,
  "DT_KCAL",  "Total energy intake, kcal/day",         "normal", NA,       0,       TRUE,
  "UPF_pct",  "Ultra-processed food intake, % energy", "normal", NA,       1,       TRUE,
  "SBP",      "Systolic blood pressure, mmHg",         "normal", NA,       1,       TRUE,
  "DBP",      "Diastolic blood pressure, mmHg",        "normal", NA,       1,       TRUE,
  "G0",       "Fasting glucose, mg/dL",                "normal", NA,       1,       TRUE,
  "USHTN",    "Hypertension, n (%)",                   "binary", "Yes",    1,       TRUE,
  "ADADM",    "Glycemic status, n (%)",                "multi",  NA,       1,       TRUE,
  "LIPIDLOW", "Lipid-lowering medicine use, n (%)",    "binary", "Yes",    1,       TRUE
)

n_b <- sum(dat$Phase == "Baseline")
n_f <- sum(dat$Phase == "Follow-up")
n_both <- sum(p4$IDNO %in% p5$IDNO)
fmt_n <- function(n) formatC(n, format = "d", big.mark = ",")

sum_normal <- function(x, d) {
  f <- paste0("%.", d, "f (%.", d, "f)")
  sprintf(f, mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}

sum_n_pct <- function(x, lvl) {
  n <- sum(x == lvl, na.rm = TRUE)
  N <- sum(!is.na(x))
  sprintf("%d (%.1f)", n, 100 * n / N)
}

p_lmer <- function(v, data) {
  dd   <- data[!is.na(data[[v]]), ]
  dd$y <- dd[[v]]
  tryCatch({
    fit <- lmerTest::lmer(y ~ Phase + (1 | FID) + (1 | ID), data = dd)
    coef(summary(fit))["PhaseFollow-up", "Pr(>|t|)"]
  }, error = function(e) {
    NA_real_
  })
}

p_gee_bin <- function(v, lvl, data) {
  dd   <- data[!is.na(data[[v]]), ] %>% arrange(FID, ID, Phase)
  dd$y <- as.integer(dd[[v]] == lvl)
  tryCatch({
    fit <- geepack::geeglm(y ~ Phase, id = FID, data = dd,
                           family = binomial, corstr = "independence")
    coef(summary(fit))["PhaseFollow-up", "Pr(>|W|)"]
  }, error = function(e) {
    NA_real_
  })
}

p_gee_multi <- function(v, data) {
  dd       <- data[!is.na(data[[v]]), ] %>% arrange(ID, Phase)
  dd$y     <- droplevels(dd[[v]])
  dd$idn   <- as.integer(dd$ID)
  dd$visit <- as.integer(dd$Phase)
  tryCatch({
    fit <- multgee::nomLORgee(y ~ Phase, data = dd, id = idn, repeated = visit,
                              LORstr = "independence")
    b   <- fit$coefficients
    idx <- grep("Phase", names(b))
    V   <- fit$robust.variance[idx, idx, drop = FALSE]
    W   <- as.numeric(t(b[idx]) %*% solve(V) %*% b[idx])
    p   <- pchisq(W, df = length(idx), lower.tail = FALSE)
    p
  }, error = function(e) {
    NA_real_
  })
}

format_p <- function(p) {
  sprintf("%.3f", p)
}

rows <- list()

for (i in seq_len(nrow(var_dict))) {
  v      <- var_dict$var[i]
  lab    <- var_dict$label[i]
  type   <- var_dict$type[i]
  d      <- var_dict$digits[i]
  do_tst <- var_dict$test[i]

  xb <- dat[[v]][dat$Phase == "Baseline"]
  xf <- dat[[v]][dat$Phase == "Follow-up"]

  if (type == "normal") {
    p <- if (do_tst) p_lmer(v, dat) else NA_real_
    rows[[length(rows) + 1]] <- data.frame(
      Characteristic = lab,
      Baseline = sum_normal(xb, d), FollowUp = sum_normal(xf, d),
      P = format_p(p))

  } else if (type == "binary") {
    lvl <- var_dict$level[i]
    p   <- if (do_tst) p_gee_bin(v, lvl, dat) else NA_real_
    rows[[length(rows) + 1]] <- data.frame(
      Characteristic = lab,
      Baseline = sum_n_pct(xb, lvl), FollowUp = sum_n_pct(xf, lvl),
      P = format_p(p))

  } else if (type == "multi") {
    p <- if (do_tst) p_gee_multi(v, dat) else NA_real_
    rows[[length(rows) + 1]] <- data.frame(
      Characteristic = lab, Baseline = "", FollowUp = "", P = format_p(p))
    for (lvl in levels(dat[[v]])) {
      rows[[length(rows) + 1]] <- data.frame(
        Characteristic = lvl,
        Baseline = sum_n_pct(xb, lvl), FollowUp = sum_n_pct(xf, lvl),
        P = "")
    }
  }
}

event_row <- function(lab, x) {
  data.frame(Characteristic = lab,
             Baseline = sprintf("%d (%.1f)", sum(x == 1, na.rm = TRUE),
                                100 * mean(x == 1, na.rm = TRUE)),
             FollowUp = "—", P = "—")
}

rows[[length(rows) + 1]] <- data.frame(
  Characteristic = "Incident events, n (%)",
  Baseline = "", FollowUp = "", P = "")
rows[[length(rows) + 1]] <- event_row("Cardiovascular disease",  p4$cvd)
rows[[length(rows) + 1]] <- event_row("Coronary heart disease",  p4$chd)
rows[[length(rows) + 1]] <- event_row("Myocardial infarction",   p4$mi)
rows[[length(rows) + 1]] <- event_row("Definite stroke",         p4$definite_stroke)

tab_df <- bind_rows(rows)

setwd(result)
write.csv(tab_df, "temp_table1_baseline_followup.csv", row.names = FALSE)

foot_1 <- paste0(
  "Continuous variables were expressed as mean (standard deviation) ",
  "and categorical variables as n (%). ",
  "Baseline and follow-up included all participants with lipidomics data at ",
  "each examination; ", fmt_n(n_both), " participants attended both examinations. ",
  "Incident events were ascertained from the baseline examination through ",
  "December 31, 2020."
)
foot_2 <- paste0(
  "P-values were obtained from linear mixed-effects models (continuous variables) ",
  "or generalized estimating equations (GEE; categorical variables), accounting for ",
  "within-family and within-person correlation."
)

ft_font <- "Times New Roman"
ft_size <- 10.5

ft <- flextable(tab_df) %>%
  set_header_labels(
    Characteristic = "",
    Baseline       = sprintf("(n=%s)", fmt_n(n_b)),
    FollowUp       = sprintf("(n=%s)", fmt_n(n_f)),
    P              = ""
  ) %>%
  add_header_row(values = c("Characteristics", "Baseline", "Follow-up", "P"),
                 colwidths = c(1, 1, 1, 1), top = TRUE) %>%
  compose(i = 1, j = "Characteristic", part = "header",
          value = as_paragraph("Characteristics", as_sup("1"))) %>%
  compose(i = 1, j = "P", part = "header",
          value = as_paragraph(as_i("P"), as_sup("2"))) %>%
  add_footer_lines(values = c(foot_1, foot_2)) %>%
  compose(i = 1, j = 1, part = "footer",
          value = as_paragraph(as_sup("1"), foot_1)) %>%
  compose(i = 2, j = 1, part = "footer",
          value = as_paragraph(as_sup("2"), foot_2)) %>%
  font(fontname = ft_font, part = "all") %>%
  fontsize(size = ft_size, part = "all") %>%
  color(color = "black", part = "all") %>%
  bg(bg = "white", part = "all") %>%
  align(align = "center", part = "header") %>%
  align(align = "center", part = "body") %>%
  align(i = 1, j = "Characteristic", align = "left", part = "header") %>%
  align(align = "left", part = "footer") %>%
  valign(valign = "top", part = "header") %>%
  line_spacing(space = 1.5, part = "body") %>%
  border_remove() %>%
  hline_top(part = "header", border = fp_border(width = 1)) %>%
  hline(i = 1, j = c("Baseline", "FollowUp"), part = "header",
        border = fp_border(width = 0.75)) %>%
  hline_bottom(part = "header", border = fp_border(width = 1)) %>%
  hline_bottom(part = "body",   border = fp_border(width = 1)) %>%
  width(j = "Characteristic", width = 2.6) %>%
  width(j = "Baseline",       width = 1.6) %>%
  width(j = "FollowUp",       width = 1.6) %>%
  width(j = "P",              width = 0.8)

setwd(tabsFigs)
title_bold  <- fp_text(bold = TRUE,  font.family = ft_font, font.size = 12)
title_plain <- fp_text(bold = FALSE, font.family = ft_font, font.size = 12)

doc <- read_docx() %>%
  body_add_fpar(fpar(
    ftext("Table 1. ", prop = title_bold),
    ftext("Characteristics of the study participants at baseline and follow-up",
          prop = title_plain)
  )) %>%
  body_add_flextable(ft)

