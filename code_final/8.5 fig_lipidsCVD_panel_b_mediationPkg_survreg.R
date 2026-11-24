# Figure 4b: mediation diagrams for lipid changes with indirect-effect P < 0.05 (main model).
# Input: result/mediation_change_mediationPkg_survreg_main.csv (5.6). Output: panel b (result/tabsFigs).

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

if (!dir.exists(tabsFigs)) dir.create(tabsFigs, recursive = TRUE)

setwd(result)
med_tag  <- "_main"
med_file <- paste0("mediation_change_mediationPkg_survreg", med_tag, ".csv")
med_all  <- read.csv(med_file, stringsAsFactors = FALSE) %>%
  mutate(exposure = "Change")

p_show  <- 0.05
med_sig <- med_all %>% filter(ACME_p < p_show)

fig_h <- 45  / 25.4

outcome_order <- c("Incident CVD", "Incident CHD", "Incident MI", "Incident definite stroke")
outcome_short <- c("Incident CVD" = "CVD", "Incident CHD" = "CHD",
                   "Incident MI"  = "MI",  "Incident definite stroke" = "Stroke")

fmt <- function(x, d) {
  x <- round(x, d); x[x == 0] <- 0
  gsub("-", "\u2212", formatC(x, format = "f", digits = d))
}
fmt_pct <- function(x) ifelse(abs(x) < 0.1 & x != 0, fmt(x, 2), fmt(x, 1))
fmt_p  <- function(p) ifelse(p < 0.001, "< 0.001", paste0("= ", formatC(p, format = "f", digits = 3)))

med_plot <- med_sig %>%
  mutate(outcome_rank = match(outcome, outcome_order),
         expo_rank    = ifelse(exposure == "Baseline", 1, 2)) %>%
  arrange(expo_rank, outcome_rank, ACME_p) %>%
  mutate(
    panel_id  = row_number(),
    out_lab   = outcome_short[outcome],
    med_lab   = ifelse(exposure == "Change", paste0("\u0394", lipids_name), lipids_name),
    head_lab  = ifelse(exposure == "Baseline", "Baseline lipid", "5-y lipid change"),
    pm_lab    = ifelse(same_sign, paste0(formatC(100 * PropMed, format = "f", digits = 1), "%\nmediated"),
                       "Inconsistent\nmediation"),
    ci_lab    = ifelse(same_sign, paste0("95% CI\n", fmt_pct(100 * PropMed_lo), "%, ",
                                         fmt_pct(100 * PropMed_hi), "%"), ""),
    pq_lab    = paste0("P ", fmt_p(ACME_p), "; FDR ", fmt_p(ACME_q))
  )
n_tri <- nrow(med_plot)

write.csv(med_plot %>% dplyr::select(panel_id, outcome, lipids_name, N_total, N_events, fu_years_median,
                                     upf_contrast, a_beta, a_lo, a_hi, b_TR, b_TR_lo, b_TR_hi,
                                     c_TR, c_TR_lo, c_TR_hi, ACME_p, ACME_q,
                                     PropMed, PropMed_lo, PropMed_hi, same_sign, covset),
          file.path(result, paste0("fig_lipidsCVD_panel_b_mediationPkg_survreg", med_tag, "_plotdata.csv")),
          row.names = FALSE)

expo_colors <- c("Baseline" = "#436891", "Change" = "#F29533")
expo_fill   <- c("Baseline" = "#E6ECF2", "Change" = "#FDEFE0")
side_fill   <- "#F0F0F0"
side_line   <- "grey45"
arrow_col   <- "grey25"

box_txt_pt  <- 7
head_pt     <- 7
in_txt_pt   <- 7
path_txt_pt <- 6.5
lab_pt      <- 7
tag_pt      <- 8

box_lwd   <- 0.35
med_lwd   <- 0.6
arrow_lwd <- 0.4

tag_mm  <- 3.5
W       <- (183 - tag_mm) / n_tri
fig_w   <- 183 / 25.4
tri_w   <- min(W, 75)
H       <- 45

box_h   <- 5.5
side_w  <- 12
x_upf   <- (W - tri_w) / 2 + 1.5 + side_w / 2
x_out   <- (W + tri_w) / 2 - 1.5 - side_w / 2
y_base  <- 14.0
x_med   <- W / 2
y_med   <- 37.3
y_head  <- 43.3
y_in    <- 25.8
y_pq    <- 7.4
y_foot  <- 2.0

text_width_mm <- function(x, size_pt, face = "plain") {
  sapply(x, function(t) {
    grid::convertWidth(
      grid::grobWidth(grid::textGrob(t, gp = grid::gpar(fontsize = size_pt,
                                                        fontfamily = "sans",
                                                        fontface = face))),
      "mm", valueOnly = TRUE)
  })
}
tmp_tiff <- tempfile(fileext = ".tiff")
if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_tiff(tmp_tiff, width = 5, height = 5, units = "in", res = 300)
} else {
  tiff(tmp_tiff, width = 5, height = 5, units = "in", res = 300)
}
med_w      <- max(text_width_mm(med_plot$med_lab, box_txt_pt, "bold")) + 4
pq_w_max   <- max(text_width_mm(med_plot$pq_lab, in_txt_pt))
foot_txt   <- paste0("")
foot_w     <- text_width_mm(foot_txt, path_txt_pt)
invisible(dev.off())
unlink(tmp_tiff)

theme_blank <- theme_void(base_family = "sans") +
  theme(plot.background  = element_rect(fill = NA, colour = NA),
        panel.background = element_rect(fill = NA, colour = NA),
        plot.margin      = margin(0, 0, 0, 0))

make_triangle <- function(r) {

  box <- function(xc, yc, w, h) c(xmin = xc - w / 2, xmax = xc + w / 2, ymin = yc - h / 2, ymax = yc + h / 2)
  b_upf <- box(x_upf, y_base, side_w, box_h)
  b_out <- box(x_out, y_base, side_w, box_h)
  b_med <- box(x_med, y_med,  med_w,  box_h)

  gap <- 0.6
  a_from <- c(x_upf + 1.5,               b_upf["ymax"] + gap)
  a_to   <- c(x_med - med_w / 2 + 2,     b_med["ymin"] - gap)
  b_from <- c(x_med + med_w / 2 - 2,     b_med["ymin"] - gap)
  b_to   <- c(x_out - 1.5,               b_out["ymax"] + gap)
  c_from <- c(b_upf["xmax"] + gap,       y_base)
  c_to   <- c(b_out["xmin"] - gap,       y_base)
  arr    <- arrow(length = unit(1.3, "mm"), type = "closed")

  ggplot() +

    annotate("segment", x = a_from[1], y = a_from[2], xend = a_to[1], yend = a_to[2],
             colour = arrow_col, linewidth = arrow_lwd, arrow = arr) +
    annotate("segment", x = b_from[1], y = b_from[2], xend = b_to[1], yend = b_to[2],
             colour = arrow_col, linewidth = arrow_lwd, arrow = arr) +
    annotate("segment", x = c_from[1], y = c_from[2], xend = c_to[1], yend = c_to[2],
             colour = arrow_col, linewidth = arrow_lwd, linetype = "dashed", arrow = arr) +

    annotate("text", x = mean(c(a_from[1], a_to[1])) - 1.6, y = mean(c(a_from[2], a_to[2])) + 0.8,
             label = "a", fontface = "italic", size = lab_pt / .pt, family = "sans", colour = "grey20") +
    annotate("text", x = mean(c(b_from[1], b_to[1])) + 1.6, y = mean(c(b_from[2], b_to[2])) + 0.8,
             label = "b", fontface = "italic", size = lab_pt / .pt, family = "sans", colour = "grey20") +
    annotate("text", x = W / 2, y = y_base - 2.1, label = "c\u2032",
             fontface = "italic", size = lab_pt / .pt, family = "sans", colour = "grey20") +

    annotate("rect", xmin = b_upf["xmin"], xmax = b_upf["xmax"], ymin = b_upf["ymin"], ymax = b_upf["ymax"],
             fill = side_fill, colour = side_line, linewidth = box_lwd) +
    annotate("rect", xmin = b_out["xmin"], xmax = b_out["xmax"], ymin = b_out["ymin"], ymax = b_out["ymax"],
             fill = side_fill, colour = side_line, linewidth = box_lwd) +
    annotate("rect", xmin = b_med["xmin"], xmax = b_med["xmax"], ymin = b_med["ymin"], ymax = b_med["ymax"],
             fill = expo_fill[[r$exposure]], colour = expo_colors[[r$exposure]], linewidth = med_lwd) +
    annotate("text", x = x_upf, y = y_base, label = "UPF", size = box_txt_pt / .pt,
             fontface = "bold", family = "sans") +
    annotate("text", x = x_out, y = y_base, label = r$out_lab, size = box_txt_pt / .pt,
             fontface = "bold", family = "sans") +
    annotate("text", x = x_med, y = y_med, label = r$med_lab, size = box_txt_pt / .pt,
             fontface = "bold", family = "sans") +

    annotate("text", x = W / 2, y = y_in + 3.0, label = r$pm_lab, size = in_txt_pt / .pt,
             fontface = "bold", lineheight = 0.85, family = "sans") +
    annotate("text", x = W / 2, y = y_in - 2.9, label = r$ci_lab, size = in_txt_pt / .pt,
             lineheight = 0.85, family = "sans", colour = "grey15") +

    annotate("text", x = W / 2, y = y_pq, label = r$pq_lab,
             size = in_txt_pt / .pt, family = "sans", colour = "grey15") +

    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, W), ylim = c(0, H), clip = "off") +
    theme_blank
}

tri_list <- lapply(seq_len(n_tri), function(i) make_triangle(med_plot[i, ]))

tag_plot <- ggplot() +
  annotate("text", x = 0, y = H, label = "b", hjust = 0, vjust = 1,
           size = tag_pt / .pt, fontface = "bold", family = "sans") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, tag_mm), ylim = c(0, H), clip = "off") +
  theme_blank

top_txt  <- "UPF \u2192 Longitudinal change in lipid \u2192 Incident CVD"
top_plot <- ggplot() +
  annotate("text", x = 1, y = y_head, label = top_txt, hjust = 0,
           size = head_pt / .pt, fontface = "bold",  family = "sans") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, W * n_tri), ylim = c(0, H), clip = "off") +
  theme_blank

foot_plot <- ggplot() +
  annotate("text", x = 1, y = y_foot, label = foot_txt, hjust = 0,
           size = path_txt_pt / .pt, fontface = "italic", colour = "grey30", family = "sans") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, W * n_tri), ylim = c(0, H), clip = "off") +
  theme_blank

design <- c(do.call(c, lapply(seq_len(n_tri + 1), function(i) patchwork::area(1, i))),
            patchwork::area(1, 2, 1, n_tri + 1),
            patchwork::area(1, 2, 1, n_tri + 1))
fig_b <- wrap_plots(c(list(tag_plot), tri_list, list(top_plot, foot_plot)), design = design,
                    widths = c(tag_mm, rep(W, n_tri)))

setwd(tabsFigs)
if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(paste0("fig_lipidsCVD_panel_b_mediationPkg_survreg", med_tag, ".tiff"), fig_b, device = ragg::agg_tiff,
         width = fig_w, height = fig_h, units = "in", dpi = 300, compression = "lzw")
} else {
  ggsave(paste0("fig_lipidsCVD_panel_b_mediationPkg_survreg", med_tag, ".tiff"), fig_b, width = fig_w, height = fig_h,
         units = "in", dpi = 300, compression = "lzw")
}
ggsave(paste0("fig_lipidsCVD_panel_b_mediationPkg_survreg", med_tag, ".pdf"), fig_b, device = cairo_pdf,
       width = fig_w, height = fig_h, units = "in")

saveRDS(list(plot = fig_b, width = fig_w, height = fig_h),
        paste0("fig_lipidsCVD_panel_b_mediationPkg_survreg", med_tag, ".rds"))
