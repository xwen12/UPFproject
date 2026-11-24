# Figure 2b-c: forest plots of total UPF and UPF subclasses with incident CVD outcomes.
# Input: 2.0.1 and 2.0.2 results. Output: panels b and c (result/tabsFigs).

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
raw_dir  <- "path/to/SHS_raw_data"
library(tidyverse)
library(ggplot2)
library(patchwork)

fig_font <- "Arial"

rawData   <- raw_dir
readyData <- file.path(proj_dir, "data")
result    <- file.path(proj_dir, "result")
tabsFigs  <- file.path(proj_dir, "result", "tabsFigs")

if (!dir.exists(tabsFigs)) dir.create(tabsFigs, recursive = TRUE)

setwd(result)
total_raw    <- read.csv("upf_total_cvd_outcomes.csv",    stringsAsFactors = FALSE)
subclass_raw <- read.csv("upf_subclass_cvd_outcomes.csv", stringsAsFactors = FALSE)

setwd(tabsFigs)
panel_a_info <- readRDS("fig1_panel_a_withChange.rds")

fig_w_total <- 183 / 25.4
fig_h_total <- 130 / 25.4
gap_in      <- 0.08

b_width  <- fig_w_total - panel_a_info$width
b_height <- panel_a_info$height
c_width  <- fig_w_total
c_height <- fig_h_total - panel_a_info$height - gap_in

outcome_short <- c("Incident CVD"             = "CVD",
                   "Incident CHD"             = "CHD",
                   "Incident MI"              = "MI",
                   "Incident definite stroke" = "Stroke")

short_labels <- c(
  "Ready-to-eat/heat mixed dishes"                   = "Ready-to-eat mixed dishes",
  "Packaged sweet snacks & desserts"                 = "Sweet snacks & desserts",
  "Sugar/artificially sweetened beverages"           = "Sweetened beverages",
  "Sauces, spreads & condiments"                     = "Sauces & spreads",
  "Meat/poultry/seafood-based ready-to-eat products" = "Ready-to-eat meat products",
  "Ultra-processed breads & breakfast foods"         = "Breads & breakfast foods",
  "Packaged savory snacks"                           = "Savory snacks",
  "Dairy-based desserts"                             = "Dairy desserts"
)

prep_forest <- function(d, x_min, x_max) {
  d %>%
    mutate(
      sig_cat = case_when(q_value < 0.05 ~ "q",
                          P_value < 0.05 ~ "p",
                          TRUE           ~ "ns"),
      lo_plot = pmax(Lower_CI, x_min),
      hi_plot = pmin(Upper_CI, x_max),
      lo_cut  = Lower_CI < x_min,
      hi_cut  = Upper_CI > x_max,
      hr_out  = HR < x_min | HR > x_max,
      ci_txt  = HR_95CI
    )
}

b_xlim   <- c(0.5, 2.25)
b_breaks <- c(0.5, 1, 2)
b_df <- total_raw %>%
  prep_forest(b_xlim[1], b_xlim[2]) %>%
  mutate(y_pos     = 5 - Outcome_order,
         name_txt  = recode(Outcome, !!!outcome_short),
         color_key = "Total UPF")

c_axis <- list(
  "Incident CVD"             = list(xlim = c(0.25,  4.5), breaks = c(0.25, 0.5, 1, 2, 4)),
  "Incident CHD"             = list(xlim = c(0.25,  4.5), breaks = c(0.25, 0.5, 1, 2, 4)),
  "Incident MI"              = list(xlim = c(0.2,   4.5), breaks = c(0.25, 0.5, 1, 2, 4)),
  "Incident definite stroke" = list(xlim = c(0.008, 4.5), breaks = c(0.01, 0.1, 1, 4))
)

c_df <- subclass_raw %>%
  mutate(x_lo = sapply(Outcome, function(o) c_axis[[o]]$xlim[1]),
         x_hi = sapply(Outcome, function(o) c_axis[[o]]$xlim[2])) %>%
  prep_forest(.$x_lo, .$x_hi) %>%
  mutate(row_id    = match(Exposure_label, names(short_labels)),
         y_pos     = 9 - row_id,
         name_txt  = recode(Exposure_label, !!!short_labels),
         color_key = Exposure_label)

nova4_colors <- c(
  "Total UPF"                                        = "#1A1A1A",
  "Packaged sweet snacks & desserts"                 = "#E0BE4F",
  "Sugar/artificially sweetened beverages"           = "#74AAA5",
  "Dairy-based desserts"                             = "#AAAAAA",
  "Ultra-processed breads & breakfast foods"         = "#F29533",
  "Sauces, spreads & condiments"                     = "#DD5454",
  "Packaged savory snacks"                           = "#436891",
  "Meat/poultry/seafood-based ready-to-eat products" = "#5FA64C",
  "Ready-to-eat/heat mixed dishes"                   = "#A7799A"
)

name_txt_pt   <- 7
ci_txt_pt     <- 7
header_pt     <- 7
subheader_pt  <- 7
axis_txt_pt   <- 7
axis_title_pt <- 7
tag_pt        <- 8

ci_lwd      <- 0.6
ref_lwd     <- 0.3
axis_lwd    <- 0.3
pt_size     <- 2.1
dia_size    <- 2.6
pt_stroke   <- 0.6
key_h_in    <- 0.18
key_col     <- "#1A1A1A"

key_df <- tibble(
  sig_cat = c("q", "p", "ns"),
  lab     = c("FDR < 0.05", "P < 0.05", "P \u2265 0.05")
)

ci_pad_in   <- 0.06

text_width_in <- function(x, size_pt, face = "plain") {
  sapply(x, function(t) {
    grid::convertWidth(
      grid::grobWidth(grid::textGrob(t, gp = grid::gpar(fontsize = size_pt,
                                                        fontfamily = fig_font,
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
w_names_c <- 1.03 * max(text_width_in(short_labels, name_txt_pt))
w_names_b <- 1.03 * max(text_width_in(outcome_short, name_txt_pt, face = "bold"))
w_ci_b    <- 1.03 * max(text_width_in(c(b_df$ci_txt, "HR (95% CI)"), ci_txt_pt))
w_xtitle  <- 1.03 * text_width_in("HR per 10% of energy intake", axis_title_pt)
key_df$w  <- 1.03 * text_width_in(key_df$lab, ci_txt_pt)
invisible(dev.off())
unlink(tmp_tiff)

outer_in   <- 0.05
tag_in     <- 0.06
ci_w_b     <- ci_pad_in + w_ci_b + 0.03
b_name_in  <- w_names_b + 0.06
b_forest   <- b_width - tag_in - outer_in - b_name_in - ci_w_b
c_name_in  <- w_names_c + 0.10
c_gap_in   <- 0.15
c_forest   <- (c_width - tag_in - outer_in - c_name_in - 3 * c_gap_in) / 4

theme_forest <- theme_classic(base_size = 7, base_family = fig_font) +
  theme(
    panel.background  = element_rect(fill = "white", colour = NA),
    panel.grid        = element_blank(),
    axis.line.x       = element_line(colour = "grey80", linewidth = axis_lwd),
    axis.line.y       = element_blank(),
    axis.ticks.y      = element_blank(),
    axis.ticks.x      = element_line(colour = "grey80", linewidth = axis_lwd),
    axis.ticks.length = unit(0.06, "cm"),
    axis.text.y       = element_blank(),
    axis.title.y      = element_blank(),
    axis.text.x       = element_text(size = axis_txt_pt, colour = "black"),
    axis.title.x      = element_text(size = axis_title_pt, colour = "black",
                                     lineheight = 0.9, margin = margin(t = 2)),
    plot.tag          = element_text(size = tag_pt, face = "bold", family = fig_font),
    plot.tag.position = "topleft",
    plot.margin       = margin(t = 2, r = 0, b = 2, l = 0)
  )

theme_textcol <- theme_forest +
  theme(
    axis.line.x  = element_blank(),
    axis.ticks.x = element_line(colour = NA),
    axis.text.x  = element_text(size = axis_txt_pt, colour = NA)
  )

make_forest <- function(d, header, x_title, xlim, breaks, y_lim, y_head, sep_y = NULL) {
  ggplot(d, aes(y = y_pos)) +

    { if (!is.null(sep_y)) geom_hline(yintercept = sep_y, colour = "grey85", linewidth = 0.25) } +

    annotate("segment", x = 1, xend = 1, y = y_lim[1], yend = y_head - 0.4,
             linetype = "dashed", colour = "grey50", linewidth = ref_lwd) +

    annotate("text", x = 1, y = y_head, label = header,
             size = header_pt / .pt, fontface = "bold", colour = "black", family = fig_font) +

    geom_segment(aes(x = lo_plot, xend = hi_plot, yend = y_pos, colour = color_key),
                 linewidth = ci_lwd) +

    geom_segment(data = d %>% filter(hi_cut),
                 aes(x = xlim[2] / 1.25, xend = xlim[2], yend = y_pos, colour = color_key),
                 arrow = arrow(length = unit(0.07, "cm"), type = "closed"), linewidth = ci_lwd) +
    geom_segment(data = d %>% filter(lo_cut),
                 aes(x = xlim[1] * 1.25, xend = xlim[1], yend = y_pos, colour = color_key),
                 arrow = arrow(length = unit(0.07, "cm"), type = "closed"), linewidth = ci_lwd) +

    geom_text(data = d %>% filter(HR < xlim[1]),
              aes(x = xlim[1] * 1.08, y = y_pos + 0.42, label = sprintf("HR %.2f", HR)),
              hjust = 0, size = 6 / .pt, colour = "grey30", family = fig_font) +
    geom_text(data = d %>% filter(HR > xlim[2]),
              aes(x = xlim[2] / 1.08, y = y_pos + 0.42, label = sprintf("HR %.2f", HR)),
              hjust = 1, size = 6 / .pt, colour = "grey30", family = fig_font) +

    geom_point(data = d %>% filter(!hr_out),
               aes(x = HR, colour = color_key,
                   fill  = ifelse(sig_cat == "ns", "hollow", color_key),
                   shape = sig_cat, size = sig_cat),
               stroke = pt_stroke) +

    scale_colour_manual(values = nova4_colors, guide = "none") +
    scale_fill_manual(values = c(nova4_colors, "hollow" = "white"), guide = "none") +
    scale_shape_manual(values = c("q" = 23, "p" = 21, "ns" = 21), guide = "none") +
    scale_size_manual(values = c("q" = dia_size, "p" = pt_size, "ns" = pt_size), guide = "none") +
    scale_x_log10(breaks = breaks, labels = format(breaks, drop0trailing = TRUE, trim = TRUE),
                  expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = xlim, ylim = y_lim, clip = "off") +
    labs(x = x_title) +
    theme_forest
}

make_ci_text <- function(d, y_lim, y_head, x_title_blank, col_w, header = "HR (95% CI)") {
  ggplot(d, aes(y = y_pos)) +
    geom_text(aes(x = ci_pad_in, label = ci_txt), hjust = 0, lineheight = 0.9,
              size = ci_txt_pt / .pt, colour = "grey15", family = fig_font) +
    annotate("text", x = ci_pad_in, y = y_head, label = header, hjust = 0, lineheight = 0.9,
             size = subheader_pt / .pt, colour = "grey30", family = fig_font) +
    scale_x_continuous(limits = c(0, col_w), expand = c(0, 0)) +
    scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(x = x_title_blank) +
    theme_textcol
}

make_key <- function(width_in, x0 = NULL) {
  key_total_in <- sum(0.07 + key_df$w + 0.16) - 0.16
  if (is.null(x0)) x0 <- (width_in - key_total_in) / 2
  k <- key_df %>%
    mutate(sym_x = x0 + c(0, cumsum(0.07 + w + 0.16)[-n()]),
           lab_x = sym_x + 0.07,
           fill  = ifelse(sig_cat == "ns", "white", key_col),
           shape = c(23, 21, 21),
           size  = c(dia_size, pt_size, pt_size))
  ggplot(k) +
    geom_point(aes(x = sym_x, y = 0.5, shape = I(shape), fill = I(fill), size = I(size)),
               colour = key_col, stroke = pt_stroke) +
    geom_text(aes(x = lab_x, y = 0.5, label = lab), hjust = 0,
              size = ci_txt_pt / .pt, colour = "grey15", family = fig_font) +
    scale_x_continuous(limits = c(0, width_in), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void(base_family = fig_font) +
    theme(plot.margin = margin(0, 0, 0, 0))
}

b_ylim <- c(0.4, 5.0)
b_head <- 4.75

b_names <- ggplot(b_df, aes(y = y_pos)) +
  geom_text(aes(x = 0, label = name_txt), hjust = 0, fontface = "bold",
            size = name_txt_pt / .pt, colour = "grey15", family = fig_font) +
  scale_x_continuous(limits = c(0, b_name_in), expand = c(0, 0)) +
  scale_y_continuous(limits = b_ylim, expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = " ") +
  theme_textcol

b_forest_plot <- make_forest(b_df, header = "Total UPF",
                             x_title = "HR per 10% of energy intake",
                             xlim = b_xlim, breaks = b_breaks,
                             y_lim = b_ylim, y_head = b_head)
b_ci_plot     <- make_ci_text(b_df, y_lim = b_ylim, y_head = b_head, x_title_blank = " ",
                              col_w = ci_w_b)

b_row <- b_names + b_forest_plot + b_ci_plot +
  plot_layout(widths = c(b_name_in, b_forest, ci_w_b))
fig1b <- b_row / make_key(b_width, x0 = tag_in + b_name_in + 0.05) +
  plot_layout(heights = c(b_height - key_h_in, key_h_in))

c_ylim <- c(0.4, 9.0)
c_head <- 8.75

if (c_forest >= w_xtitle + 0.05) {
  c_xtitle       <- "HR per 10% of energy intake"
  c_xtitle_blank <- " "
} else {
  c_xtitle       <- "HR per 10% of\nenergy intake"
  c_xtitle_blank <- " \n "
}

c_name_df <- c_df %>% distinct(name_txt, y_pos)
c_names <- ggplot(c_name_df, aes(y = y_pos)) +
  geom_text(aes(x = c_name_in - 0.06, label = name_txt), hjust = 1,
            size = name_txt_pt / .pt, colour = "grey15", family = fig_font) +
  scale_x_continuous(limits = c(0, c_name_in), expand = c(0, 0)) +
  scale_y_continuous(limits = c_ylim, expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = c_xtitle_blank) +
  theme_textcol

c_plots <- list(c_names)
for (k in seq_along(outcome_short)) {
  oc <- names(outcome_short)[k]
  d  <- c_df %>% filter(Outcome == oc)
  if (k > 1) c_plots <- c(c_plots, list(plot_spacer()))
  c_plots <- c(c_plots, list(
    make_forest(d, header = outcome_short[[oc]], x_title = c_xtitle,
                xlim = c_axis[[oc]]$xlim, breaks = c_axis[[oc]]$breaks,
                y_lim = c_ylim, y_head = c_head)
  ))
}

c_row <- wrap_plots(c_plots, nrow = 1) +
  plot_layout(widths = c(c_name_in, c_forest, rep(c(c_gap_in, c_forest), 3)))
fig1c <- c_row / make_key(c_width) +
  plot_layout(heights = c(c_height - key_h_in, key_h_in))

setwd(tabsFigs)
ggsave(fig1b, file = "fig1_panel_b_totalUPF.tiff", bg = "white",
       width = b_width, height = b_height, units = "in", dpi = 600, compression = "lzw",
       device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_tiff else "tiff")
ggsave(fig1c, file = "fig1_panel_c_subclass.tiff", bg = "white",
       width = c_width, height = c_height, units = "in", dpi = 600, compression = "lzw",
       device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_tiff else "tiff")

saveRDS(list(plot = patchwork::patchworkGrob(fig1b), width = b_width, height = b_height),
        file = "fig1_panel_b_totalUPF.rds")
saveRDS(list(plot = patchwork::patchworkGrob(fig1c), width = c_width, height = c_height),
        file = "fig1_panel_c_subclass.rds")
