# Figure 4a: forest plots of baseline lipid levels and 5-year lipid changes with incident outcomes.
# Input: 3.1-3.3 and 4.1-4.3 results. Output: panel a (result/tabsFigs).

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
raw_dir  <- "path/to/SHS_raw_data"
library(tidyverse)
library(ggplot2)
library(patchwork)

rawData   <- raw_dir
readyData <- file.path(proj_dir, "data")
result    <- file.path(proj_dir, "result")
tabsFigs  <- file.path(proj_dir, "result", "tabsFigs")

lipid_info_path <- file.path(raw_dir, "lipidomics data", "lipid_4_6_compare.csv")

if (!dir.exists(tabsFigs)) dir.create(tabsFigs, recursive = TRUE)

exposure_files <- c("longitu" = "Baseline", "change" = "Change")
outcome_files  <- c("cvd" = "CVD", "chd" = "CHD", "mi" = "MI", "stroke" = "Stroke")
set_files      <- c("riskLipids" = "UPF-positive", "protLipids" = "UPF-inverse")

setwd(result)
res_long <- expand_grid(ex = names(exposure_files),
                        out = names(outcome_files),
                        st = names(set_files)) %>%
  pmap_dfr(function(ex, out, st) {
    f <- paste0(ex, "_", out, "_", st, ".csv")
    read.csv(f, stringsAsFactors = FALSE) %>%
      dplyr::select(lipids_name, beta, se, HR, lower_CI, upper_CI, pvals, qvals) %>%
      mutate(exposure  = exposure_files[[ex]],
             outcome   = outcome_files[[out]],
             lipid_set = set_files[[st]],
             file      = f)
  })

res_long_out <- res_long %>%
  mutate(HR_95CI = sprintf("%.2f (%.2f, %.2f)", HR, lower_CI, upper_CI)) %>%
  arrange(lipid_set, lipids_name, exposure, outcome)
write.csv(res_long_out, file.path(result, "lipids42_cvdOutcomes_bsl_chg_long.csv"), row.names = FALSE)

fig_w <- 183 / 25.4
fig_h <- 90  / 25.4

sel_lipids <- res_long %>%
  group_by(lipids_name, lipid_set) %>%
  summarise(n_p05   = sum(pvals < 0.05),
            n_q05   = sum(qvals < 0.05),
            min_p   = min(pvals),
            .groups = "drop") %>%
  filter(n_p05 >= 1)

lipid_info <- read.csv(lipid_info_path, header = TRUE, stringsAsFactors = FALSE,
                       check.names = FALSE) %>%
  dplyr::select(lipids_name, Class_short) %>%
  distinct()
lipid_info[lipid_info$Class_short == "N/A",                       "Class_short"] <- "Unknown"
lipid_info[lipid_info$Class_short %in% c("FAHFA", "scFA", "USFA"), "Class_short"] <- "FA"
lipid_info[lipid_info$Class_short %in% c("Glc"),                  "Class_short"] <- "GlcCer"

cls_order <- c("AC", "CE", "CER", "DAG", "FA", "GlcCer", "PC", "PE", "PG", "PI",
               "PS", "SM", "TAG", "Unknown")

row_info <- sel_lipids %>%
  left_join(lipid_info, by = "lipids_name") %>%
  rename(lipid_class = Class_short) %>%
  mutate(class_rank = match(lipid_class, cls_order),
         label      = lipids_name) %>%
  arrange(lipid_set, class_rank, lipids_name) %>%
  group_by(lipid_set) %>%
  mutate(row_id = row_number(),
         n_rows = n(),
         y_pos  = n_rows - row_id + 1) %>%
  ungroup()

x_lim    <- c(0.2, 3.3)
x_breaks <- c(0.25, 0.5, 1, 2)

dodge <- 0.24

plot_df <- res_long %>%
  inner_join(row_info %>% dplyr::select(lipids_name, lipid_set, label, lipid_class, y_pos),
             by = c("lipids_name", "lipid_set")) %>%
  mutate(
    exposure = factor(exposure, levels = c("Baseline", "Change")),
    outcome  = factor(outcome,  levels = outcome_files),
    y_plot   = y_pos + ifelse(exposure == "Baseline", dodge, -dodge),
    sig_cat  = case_when(qvals < 0.05 ~ "q",
                         pvals < 0.05 ~ "p",
                         TRUE         ~ "ns"),
    fill_key = ifelse(sig_cat == "ns", "hollow", as.character(exposure)),
    lo_plot  = pmax(lower_CI, x_lim[1]),
    hi_plot  = pmin(upper_CI, x_lim[2]),
    lo_cut   = lower_CI < x_lim[1],
    hi_cut   = upper_CI > x_lim[2],
    hr_out   = HR < x_lim[1] | HR > x_lim[2]
  )

write.csv(plot_df %>% dplyr::select(lipid_set, lipids_name, label, exposure, outcome,
                             HR, lower_CI, upper_CI, pvals, qvals, sig_cat) %>%
            arrange(lipid_set, label, exposure, outcome),
          file.path(result, "fig_lipidsCVD_panel_a_plotdata.csv"), row.names = FALSE)

expo_colors <- c("Baseline" = "#436891", "Change" = "#F29533")
expo_labels <- c("Baseline" = "Baseline lipid level",
                 "Change"   = "Longitudinal change in lipid")

band_fill <- "#F0F0F0"

name_txt_pt   <- 7
header_pt     <- 7
title_pt      <- 7
axis_txt_pt   <- 7
axis_title_pt <- 7
key_txt_pt    <- 7
tag_pt        <- 8

ci_lwd    <- 0.4
ref_lwd   <- 0.3
axis_lwd  <- 0.3
sep_lwd   <- 0.25
pt_size   <- 1.3
dia_size  <- 1.8
pt_stroke <- 0.5
key_col   <- "#1A1A1A"

text_width_in <- function(x, size_pt, face = "plain") {
  sapply(x, function(t) {
    grid::convertWidth(
      grid::grobWidth(grid::textGrob(t, gp = grid::gpar(fontsize = size_pt,
                                                        fontfamily = "sans",
                                                        fontface = face))),
      "in", valueOnly = TRUE)
  })
}
tmp_tiff <- tempfile(fileext = ".tiff")
if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_tiff(tmp_tiff, width = 5, height = 5, units = "in", res = 300)
} else {
  tiff(tmp_tiff, width = 5, height = 5, units = "in", res = 300)
}
w_names   <- 1.03 * max(text_width_in(row_info$label, name_txt_pt))
w_tag     <- text_width_in("a", tag_pt, face = "bold")
key_lab_w <- 1.03 * text_width_in(c(expo_labels, "FDR < 0.05", "P < 0.05", "P ≥ 0.05"),
                                  key_txt_pt)
invisible(dev.off())
unlink(tmp_tiff)

tag_in    <- w_tag + 0.04
name_in   <- w_names + 0.10
gap_in    <- 0.15
outer_in  <- 0.05
forest_in <- (fig_w - tag_in - name_in - 3 * gap_in - outer_in) / 4

n_pos     <- sum(row_info$lipid_set == "UPF-positive")
n_inv     <- sum(row_info$lipid_set == "UPF-inverse")
t_head    <- 0.55
t_title1  <- 1.45
t_rows1   <- t_title1 + 0.95 + (seq_len(n_pos) - 1)
t_title2  <- max(t_rows1) + 1.15
t_rows2   <- t_title2 + 0.95 + (seq_len(n_inv) - 1)
t_axis    <- max(t_rows2) + 0.60
below_in  <- 0.46
row_in    <- (fig_h - below_in - 4 / 72) / t_axis
t_bottom  <- t_axis + below_in / row_in
y_lim     <- c(-t_bottom, 0)

axis_lab    <- format(x_breaks, drop0trailing = TRUE, trim = TRUE)
y_axis      <- -t_axis
tick_len    <- 0.06 / 2.54 / row_in
axis_num_dy <- 0.07 / row_in
axis_tit_dy <- 0.20 / row_in
y_key       <- y_axis - 0.40 / row_in

theme_blank <- theme_void(base_size = 7, base_family = "sans") +
  theme(
    panel.background = element_rect(fill = NA, colour = NA),
    plot.background  = element_rect(fill = NA, colour = NA),
    plot.margin      = margin(t = 2, r = 0, b = 2, l = 0)
  )

make_forest <- function(d, header, blk_y) {
  y_axis <- y_axis; tick_len <- tick_len; axis_num_dy <- axis_num_dy; x_lim <- x_lim
  ticks <- tibble(x = x_breaks, lab = axis_lab)
  ggplot(d, aes(y = y_plot)) +

    annotate("segment", x = x_lim[1], xend = x_lim[2], y = blk_y, yend = blk_y,
             colour = "grey75", linewidth = sep_lwd) +

    annotate("segment", x = 1, xend = 1, y = y_axis, yend = -t_head - 0.45,
             linetype = "dashed", colour = "grey50", linewidth = ref_lwd) +

    annotate("text", x = 1, y = -t_head, label = header,
             size = header_pt / .pt, fontface = "bold", colour = "black", family = "sans") +

    annotate("segment", x = x_lim[1], xend = x_lim[2], y = y_axis, yend = y_axis,
             colour = "grey80", linewidth = axis_lwd) +
    geom_segment(data = ticks, aes(x = x, xend = x, y = y_axis, yend = y_axis - tick_len),
                 inherit.aes = FALSE, colour = "grey80", linewidth = axis_lwd) +
    geom_text(data = ticks, aes(x = x, y = y_axis - axis_num_dy, label = lab),
              inherit.aes = FALSE, vjust = 1, size = axis_txt_pt / .pt,
              colour = "black", family = "sans") +
    annotate("text", x = sqrt(x_lim[1] * x_lim[2]), y = y_axis - axis_tit_dy,
             label = "HR per SD", vjust = 1, size = axis_title_pt / .pt,
             colour = "black", family = "sans") +

    geom_segment(aes(x = lo_plot, xend = hi_plot, yend = y_plot, colour = exposure),
                 linewidth = ci_lwd) +

    geom_segment(data = d %>% filter(hi_cut),
                 aes(x = x_lim[2] / 1.25, xend = x_lim[2], yend = y_plot, colour = exposure),
                 arrow = arrow(length = unit(0.06, "cm"), type = "closed"), linewidth = ci_lwd) +
    geom_segment(data = d %>% filter(lo_cut),
                 aes(x = x_lim[1] * 1.25, xend = x_lim[1], yend = y_plot, colour = exposure),
                 arrow = arrow(length = unit(0.06, "cm"), type = "closed"), linewidth = ci_lwd) +

    geom_point(data = d %>% filter(!hr_out),
               aes(x = HR, colour = exposure, fill = fill_key, shape = sig_cat, size = sig_cat),
               stroke = pt_stroke) +

    scale_colour_manual(values = expo_colors, guide = "none") +
    scale_fill_manual(values = c(expo_colors, "hollow" = "white"), guide = "none") +
    scale_shape_manual(values = c("q" = 23, "p" = 21, "ns" = 21), guide = "none") +
    scale_size_manual(values = c("q" = dia_size, "p" = pt_size, "ns" = pt_size), guide = "none") +
    scale_x_log10(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = x_lim, ylim = y_lim, clip = "off") +
    theme_blank
}

make_names <- function(r, titles) {
  ggplot(r, aes(y = y_pos)) +
    geom_text(aes(x = 0, label = label), hjust = 0,
              size = name_txt_pt / .pt, colour = "grey15", family = "sans") +
    geom_text(data = titles, aes(x = 0, y = y, label = lab), hjust = 0,
              size = title_pt / .pt, fontface = "bold", colour = "black", family = "sans") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, name_in), ylim = y_lim, clip = "off") +
    theme_blank
}

make_tag <- function(tag = "a") {
  ggplot() +
    annotate("text", x = 0, y = 0, label = tag, hjust = 0, vjust = 1,
             size = tag_pt / .pt, fontface = "bold", family = "sans") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, tag_in), ylim = y_lim, clip = "off") +
    theme_blank
}

make_key <- function(width_in) {
  y_key <- y_key
  sym_w <- 0.20; lab_gap <- 0.05; item_gap <- 0.16; group_gap <- 0.30
  k <- tibble(lab   = c(expo_labels, "FDR < 0.05", "P < 0.05", "P ≥ 0.05"),
              w     = key_lab_w,
              shape = c(21, 21, 23, 21, 21),
              fill  = c(expo_colors, key_col, key_col, "white"),
              col   = c(expo_colors, key_col, key_col, key_col),
              size  = c(pt_size, pt_size, dia_size, pt_size, pt_size),
              bar   = c(TRUE, TRUE, FALSE, FALSE, FALSE),
              extra = c(0, 0, group_gap, 0, 0))
  step <- sym_w + lab_gap + k$w + item_gap
  tot  <- sum(step) - item_gap + sum(k$extra)
  k <- k %>% mutate(sym_x = (width_in - tot) / 2 + c(0, cumsum(step)[-n()]) +
                      cumsum(extra) + sym_w / 2,
                    lab_x = sym_x + sym_w / 2 + lab_gap)
  ggplot(k) +
    geom_segment(data = k %>% filter(bar),
                 aes(x = sym_x - sym_w / 2, xend = sym_x + sym_w / 2, y = y_key, yend = y_key,
                     colour = I(col)), linewidth = ci_lwd) +
    geom_point(aes(x = sym_x, y = y_key, shape = I(shape), fill = I(fill), colour = I(col),
                   size = I(size)), stroke = pt_stroke) +
    geom_text(aes(x = lab_x, y = y_key, label = lab), hjust = 0,
              size = key_txt_pt / .pt, colour = "grey15", family = "sans") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, width_in), ylim = y_lim, clip = "off") +
    theme_blank
}

make_bands <- function(bands, width_in) {
  ggplot(bands) +
    geom_rect(aes(xmin = -0.03, xmax = width_in, ymin = ymin, ymax = ymax),
              fill = band_fill, colour = NA) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, width_in), ylim = y_lim, clip = "off") +
    theme_blank
}

row_info <- row_info %>%
  mutate(t     = ifelse(lipid_set == "UPF-positive", t_rows1[row_id], t_rows2[row_id]),
         y_pos = -t)
plot_df  <- plot_df %>%
  dplyr::select(-y_pos) %>%
  left_join(row_info %>% dplyr::select(lipids_name, lipid_set, y_pos), by = c("lipids_name", "lipid_set")) %>%
  mutate(y_plot = y_pos + ifelse(exposure == "Baseline", dodge, -dodge))

bands <- row_info %>%
  group_by(lipid_set) %>% arrange(row_id, .by_group = TRUE) %>%
  mutate(run = cumsum(lipid_class != lag(lipid_class, default = first(lipid_class)))) %>%
  group_by(lipid_set, run, lipid_class) %>%
  summarise(ymin = min(y_pos) - 0.5, ymax = max(y_pos) + 0.5, .groups = "drop") %>%
  mutate(shaded = run %% 2 == 0)
bands <- bands %>% filter(shaded)
blk_y <- -(t_title2 - 0.55)

titles <- tibble(y   = c(-t_title1, -t_title2),
                 lab = c("Lipids positively associated with UPF",
                         "Lipids inversely associated with UPF"))

forests  <- map(levels(plot_df$outcome), function(o)
  make_forest(plot_df %>% filter(outcome == o), header = o, blk_y = blk_y))
band_w     <- name_in + 4 * forest_in + 3 * gap_in
band_plot  <- make_bands(bands, band_w)
key_plot   <- make_key(band_w)

widths <- c(tag_in, name_in, forest_in, gap_in, forest_in, gap_in, forest_in, gap_in, forest_in,
            outer_in)
gap <- function() ggplot() + coord_cartesian(ylim = y_lim) + theme_blank
cells  <- list(make_tag("a"), make_names(row_info, titles),
               forests[[1]], gap(), forests[[2]], gap(),
               forests[[3]], gap(), forests[[4]], gap())
design <- c(
  patchwork::area(1, 2, 1, 9),
  do.call(c, lapply(seq_along(widths), function(i) patchwork::area(1, i))),
  patchwork::area(1, 2, 1, 9)
)
fig_a <- wrap_plots(c(list(band_plot), cells, list(key_plot)), design = design, widths = widths)

setwd(tabsFigs)
if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave("fig_lipidsCVD_panel_a_forest.tiff", fig_a, device = ragg::agg_tiff,
         width = fig_w, height = fig_h, units = "in", dpi = 300, compression = "lzw")
} else {
  ggsave("fig_lipidsCVD_panel_a_forest.tiff", fig_a, width = fig_w, height = fig_h,
         units = "in", dpi = 300, compression = "lzw")
}
ggsave("fig_lipidsCVD_panel_a_forest.pdf", fig_a, device = cairo_pdf,
       width = fig_w, height = fig_h, units = "in")

saveRDS(list(plot = fig_a, width = fig_w, height = fig_h),
        "fig_lipidsCVD_panel_a_forest.rds")
