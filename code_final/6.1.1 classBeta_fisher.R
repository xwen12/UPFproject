# Figure 3b: total UPF betas by lipid class, with over-representation tests by class.
# Input: result/P4P5all.csv (2.1), lipid class annotation. Output: panel b (result/tabsFigs).

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
raw_dir  <- "path/to/SHS_raw_data"

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

fig_font <- "Arial"

result   <- file.path(proj_dir, "result")
tabsFigs <- file.path(proj_dir, "result", "tabsFigs")
lipid_info_path <- file.path(raw_dir, "lipidomics data", "lipid_4_6_compare.csv")

min_n <- 5

setwd(result)
P4P5all <- read.table("P4P5all.csv", sep = ",", header = TRUE, check.names = FALSE)

lipid_info <- read.csv(lipid_info_path, header = TRUE, stringsAsFactors = FALSE,
                       check.names = FALSE)

dt0 <- P4P5all %>%
  filter(known) %>%
  select(lipids_name, beta, se, pvals, qvals) %>%
  left_join(lipid_info[, c("lipids_name", "Class_short")], by = "lipids_name")

dt0[dt0$Class_short == "N/A",                       "Class_short"] <- "Unknown"
dt0[dt0$Class_short %in% c("FAHFA", "scFA", "USFA"), "Class_short"] <- "FA"
dt0[dt0$Class_short %in% c("Glc"),                  "Class_short"] <- "GlcCer"
dt0 <- dt0 %>% filter(Class_short != "Unknown")

dt0 <- dt0 %>%
  mutate(sig_grp = case_when(qvals < 0.05 & beta > 0 ~ "up",
                             qvals < 0.05 & beta < 0 ~ "down",
                             TRUE                    ~ "ns"))

cls_n <- dt0 %>% count(Class_short, name = "n_lipids")

cls_order <- sort(unique(cls_n$Class_short), method = "radix")

cls_colors <- c(
  AC = "#3F7F33", CE = "#2F4F73", CER = "#C94C7C", DAG = "#4C78A8",
  FA = "#B8860B", GlcCer = "#4D4D4D", PC = "#2C6DA4", PE = "#6A5AAE",
  PG = "#C75D1A", PI = "#AD3F7A", PS = "#8E6BBF", SM = "#1F7A6E",
  TAG = "#A0522D")
miss_col <- setdiff(cls_order, names(cls_colors))
cls_colors <- c(cls_colors, setNames(rep("#7A7A7A", length(miss_col)), miss_col))

test_list <- list()
for (cls in cls_order) {
  in_cls <- dt0$Class_short == cls
  sig    <- dt0$qvals < 0.05

  tab <- matrix(c(sum(in_cls & sig),  sum(in_cls & !sig),
                  sum(!in_cls & sig), sum(!in_cls & !sig)),
                nrow = 2, byrow = TRUE)
  tested <- sum(in_cls) >= min_n

  expected <- outer(rowSums(tab), colSums(tab)) / sum(tab)
  min_exp  <- min(expected)

  if (tested) {
    if (min_exp >= 5) {
      pt        <- prop.test(x = tab[, 1], n = rowSums(tab),
                             alternative = "greater", correct = FALSE)
      p_enrich  <- pt$p.value
      test_used <- "chi-square"
    } else {
      p_enrich  <- fisher.test(tab, alternative = "greater")$p.value
      test_used <- "Fisher"
    }
    ft2 <- fisher.test(tab, alternative = "two.sided")
  }

  test_list[[cls]] <- data.frame(
    Class_short   = cls,
    n_lipids      = sum(in_cls),
    n_up          = sum(in_cls & dt0$sig_grp == "up"),
    n_down        = sum(in_cls & dt0$sig_grp == "down"),
    pct_sig       = 100 * mean(sig[in_cls]),
    median_beta   = median(dt0$beta[in_cls]),
    tested        = tested,
    min_expected  = min_exp,
    test_used     = if (tested) test_used          else NA_character_,
    OR            = if (tested) unname(ft2$estimate) else NA_real_,
    OR_lo         = if (tested) ft2$conf.int[1]      else NA_real_,
    OR_hi         = if (tested) ft2$conf.int[2]      else NA_real_,
    enrich_p      = if (tested) p_enrich           else NA_real_)
}
class_res <- bind_rows(test_list)

write.table(class_res, "classBeta_fisher.csv", sep = ",",
            row.names = F, col.names = T)

axis_txt_pt  <- 7
axis_tit_pt  <- 7
tbl_txt_pt   <- 7
tag_pt       <- 8

col_pos <- "#C0392B"
col_neg <- "#436891"
col_ns  <- "grey55"

fmt_p <- function(p) {
  txt <- ifelse(p < 0.001, "<0.001",
                ifelse(p < 0.01, sprintf("%.3f", p), sprintf("%.2f", p)))
  ifelse(is.na(p), "\u2013", txt)
}

plot_df <- dt0 %>%
  filter(Class_short %in% cls_order) %>%
  mutate(Class_short = factor(Class_short, levels = rev(cls_order)),
         sig_grp     = factor(sig_grp, levels = c("ns", "down", "up"))) %>%
  arrange(sig_grp)

class_res <- class_res %>%
  mutate(Class_short = factor(Class_short, levels = rev(cls_order)),
         cls_txt     = as.character(Class_short),
         fis_txt     = fmt_p(enrich_p))

y_cls <- scale_y_discrete(drop = FALSE, expand = c(0, 0.5))

if (requireNamespace("ggbeeswarm", quietly = TRUE)) {
  dot_layer <- ggbeeswarm::geom_quasirandom(
    aes(colour = sig_grp, size = sig_grp, alpha = sig_grp), shape = 16,
    width = 0.38, orientation = "y")
} else {
  dot_layer <- geom_jitter(aes(colour = sig_grp, size = sig_grp, alpha = sig_grp),
                           shape = 16, width = 0, height = 0.3, seed = 42)
}

p_beta <- ggplot(plot_df, aes(x = beta, y = Class_short)) +

  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3, linetype = "dashed") +

  geom_boxplot(outlier.shape = NA, fill = NA, colour = "grey35",
               width = 0.6, linewidth = 0.3, fatten = 1.5, orientation = "y") +

  dot_layer +
  scale_colour_manual(values = c(up = col_pos, down = col_neg, ns = col_ns),
                    breaks = c("down", "up", "ns"),
                    labels = c("FDR < 0.05, \u03B2 < 0", "FDR < 0.05, \u03B2 > 0",
                               "FDR \u2265 0.05"),
                    name = NULL) +
  scale_size_manual(values = c(up = 1.5, down = 1.5, ns = 0.9), guide = "none") +
  scale_alpha_manual(values = c(up = 0.8, down = 0.8, ns = 0.55), guide = "none") +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 1.8, alpha = 0.8))) +

  y_cls +
  scale_x_continuous(labels = function(x) ifelse(abs(x) < 1e-9, "0",
                                                 gsub("-", "\u2212", sprintf("%.2f", x)))) +
  labs(x = "Standardized \u03B2 of total UPF", y = NULL) +

  theme_classic(base_size = 7, base_family = fig_font) +
  theme(
    axis.line.x        = element_line(colour = "grey80", linewidth = 0.3),
    axis.line.y        = element_blank(),
    axis.ticks.x       = element_line(colour = "grey80", linewidth = 0.3),
    axis.ticks.y       = element_blank(),
    axis.text.y        = element_blank(),
    axis.text.x        = element_text(size = axis_txt_pt, colour = "black"),
    axis.title.x       = element_text(size = axis_tit_pt, colour = "black",
                                      margin = margin(t = 3)),
    panel.grid.major.y = element_line(colour = "grey94", linewidth = 0.25),
    legend.position    = "bottom",
    legend.justification = "left",
    legend.text        = element_text(size = axis_txt_pt, colour = "black"),
    legend.key.size    = unit(0.3, "cm"),
    legend.margin      = margin(0, 0, 0, 0),
    legend.box.spacing = unit(2, "pt"),
    plot.margin        = margin(t = 2, r = 2, b = 2, l = 0)
  )

col_theme <- theme_void(base_size = 7, base_family = fig_font) +
  theme(plot.margin = margin(t = 0, r = 1, b = 2, l = 1))

make_col <- function(lab_col, txt_col = "black") {
  ggplot(class_res, aes(x = 0.5, y = Class_short)) +
    geom_text(aes(label = .data[[lab_col]]), colour = txt_col,
              size = tbl_txt_pt / .pt, family = fig_font) +
    y_cls +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    col_theme
}

c_class <- ggplot(class_res, aes(x = 0.5, y = Class_short)) +
  geom_tile(aes(fill = Class_short), width = 1, height = 0.86) +
  scale_fill_manual(values = cls_colors, guide = "none") +
  geom_text(aes(label = cls_txt), colour = "white",
            size = 6 / .pt, family = fig_font) +
  y_cls +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  col_theme

c_n   <- make_col("n_lipids")
c_fis <- make_col("fis_txt")

cnt_df <- bind_rows(
  class_res %>% transmute(Class_short, x = 1, lab = n_down, grp = "down"),
  class_res %>% transmute(Class_short, x = 2, lab = n_up,   grp = "up"))

c_updown <- ggplot(cnt_df, aes(x = x, y = Class_short)) +
  geom_text(aes(label = lab, colour = grp), size = tbl_txt_pt / .pt, family = fig_font) +
  scale_colour_manual(values = c(up = col_pos, down = col_neg), guide = "none") +
  y_cls +
  scale_x_continuous(limits = c(0.5, 2.5), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  col_theme

head_fill <- "grey92"
head_theme <- theme_void(base_size = 7, base_family = fig_font) +
  theme(plot.margin = margin(t = 2, r = 1, b = 1, l = 1))

make_head <- function(txt, shade = TRUE) {
  g <- ggplot() + xlim(0, 1) + ylim(0, 1)
  if (shade) g <- g + annotate("rect", xmin = 0.03, xmax = 0.97, ymin = 0, ymax = 1,
                               fill = head_fill, colour = NA)
  g + annotate("text", x = 0.5, y = 0.5, label = txt, lineheight = 0.9,
               size = tbl_txt_pt / .pt, family = fig_font) +
    coord_cartesian(expand = FALSE) + head_theme
}

h_class <- make_head("Lipid\nclass", shade = FALSE)
h_beta  <- ggplot() + theme_void()
h_n     <- make_head("Lipids,\nn")
h_p     <- make_head("P for over-\nrepresentation")

h_updown <- ggplot() + xlim(0.5, 2.5) + ylim(0, 1) +
  annotate("rect", xmin = 0.53, xmax = 2.47, ymin = 0, ymax = 1,
           fill = head_fill, colour = NA) +
  annotate("text", x = 1.5, y = 0.73, label = "FDR < 0.05, n",
           size = tbl_txt_pt / .pt, family = fig_font) +
  annotate("text", x = c(1, 2), y = 0.27, label = c("\u03B2 < 0", "\u03B2 > 0"),
           size = tbl_txt_pt / .pt, family = fig_font) +
  coord_cartesian(expand = FALSE) + head_theme

col_w_mm <- c(class = 14, beta = 25, n = 9, updown = 19, p = 19)
head_h   <- 0.075
fig_classBeta <- wrap_plots(h_class, h_beta, h_n, h_updown, h_p,
                            c_class, p_beta, c_n, c_updown, c_fis,
                            nrow = 2, widths = col_w_mm, heights = c(head_h, 1),
                            guides = "collect") &
  theme(legend.position = "bottom", legend.justification = "center")

fig_width  <- 84 / 25.4
fig_height <- 120 / 25.4

setwd(tabsFigs)
ggsave(fig_classBeta, file = "fig_classBeta_lipids.tiff",
       width = fig_width, height = fig_height, units = "in", dpi = 600,
       compression = "lzw", bg = "white",
       device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_tiff else "tiff")
ggsave(fig_classBeta, file = "fig_classBeta_lipids.pdf",
       width = fig_width, height = fig_height, units = "in", device = cairo_pdf)
saveRDS(list(plot = fig_classBeta, width = fig_width, height = fig_height),
        "fig_classBeta_lipids.rds")

