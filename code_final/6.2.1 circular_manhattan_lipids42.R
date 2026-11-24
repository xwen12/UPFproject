# Figure 3a: circular Manhattan plot of the 42 UPF-related lipids.
# Input: 2.1 results, lipid class annotation. Output: panel a (result/tabsFigs).

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
raw_dir  <- "path/to/SHS_raw_data"

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(scales)

fig_font <- "Arial"

result   <- file.path(proj_dir, "result")
tabsFigs <- file.path(proj_dir, "result", "tabsFigs")
lipid_info_path <- file.path(raw_dir, "lipidomics data", "lipid_4_6_compare.csv")

setwd(result)

upf_lipids <- read.table("UPF_lipids_all.csv", sep = ",", header = TRUE,
                         check.names = FALSE)

P4P5all  <- read.table("P4P5all.csv", sep = ",", header = TRUE, check.names = FALSE)
fdr_pval <- max(P4P5all$pvals[P4P5all$qvals < 0.05])
fdr_logP <- -log10(fdr_pval)

lipid_info <- read.csv(lipid_info_path, header = TRUE, stringsAsFactors = FALSE,
                       check.names = FALSE)

dt0 <- upf_lipids %>%
  select(lipids_name, beta, se, lower_CI, upper_CI, pvals, qvals) %>%
  left_join(lipid_info[, c("lipids_name", "Class_short")], by = "lipids_name")

dt0[dt0$Class_short == "N/A",                       "Class_short"] <- "Unknown"
dt0[dt0$Class_short %in% c("FAHFA", "scFA", "USFA"), "Class_short"] <- "FA"
dt0[dt0$Class_short %in% c("Glc"),                  "Class_short"] <- "GlcCer"

cls_order <- c("AC", "CE", "CER", "DAG", "FA", "GlcCer", "PC", "PE", "PG", "PI",
               "PS", "SM", "TAG", "Unknown")

cls_colors <- c(
  AC = "#3F7F33", CE = "#2F4F73", CER = "#C94C7C", DAG = "#4C78A8",
  FA = "#B8860B", GlcCer = "#4D4D4D", PC = "#2C6DA4", PE = "#6A5AAE",
  PG = "#C75D1A", PI = "#AD3F7A", PS = "#8E6BBF", SM = "#1F7A6E",
  TAG = "#A0522D", Unknown = "#7A7A7A")

cls_full <- c(
  AC = "AC: Acylcarnitine", CE = "CE: Cholesteryl ester", CER = "CER: Ceramide",
  DAG = "DAG: Diacylglycerol", FA = "FA: Fatty acid",
  GlcCer = "GlcCer: Glucosylceramide", PC = "PC: Phosphatidylcholine",
  PE = "PE: Phosphatidylethanolamine", PG = "PG: Phosphatidylglycerol",
  PI = "PI: Phosphatidylinositol", PS = "PS: Phosphatidylserine",
  SM = "SM: Sphingomyelin", TAG = "TAG: Triacylglycerol", Unknown = "Unknown")

cls_present <- cls_order[cls_order %in% dt0$Class_short]

plot_data <- dt0 %>%
  mutate(logP        = -log10(pvals),
         Class_short = factor(Class_short, levels = cls_present)) %>%
  arrange(Class_short, lipids_name) %>%
  mutate(id = row_number())

n_lipids    <- nrow(plot_data)
gap_radians <- 0.75
usable_arc  <- 2 * pi - gap_radians
bar_width   <- usable_arc / n_lipids
half_w      <- bar_width / 2

plot_data <- plot_data %>%
  mutate(angle       = pi / 2 - gap_radians / 2 - bar_width / 2 - (id - 1) * bar_width,
         hjust       = ifelse(cos(angle) >= 0, 0, 1),
         label_angle = ifelse(cos(angle) >= 0, angle * 180 / pi, angle * 180 / pi + 180))

inner_r  <- 5
max_bar  <- 8
max_logP <- ceiling(max(plot_data$logP, na.rm = TRUE))
to_r     <- function(lp) inner_r + (lp / max_logP) * max_bar

plot_data <- plot_data %>% mutate(bar_len = to_r(logP))

char_w <- 0.30
label_ends <- plot_data %>%
  mutate(r_end = bar_len + 0.3 + nchar(lipids_name) * char_w,
         x = r_end * cos(angle), y = r_end * sin(angle))

make_bar_poly <- function(angle, inner, outer, half_angle) {
  a <- seq(angle - half_angle, angle + half_angle, length.out = 6)
  data.frame(x = c(inner * cos(a), rev(outer * cos(a))),
             y = c(inner * sin(a), rev(outer * sin(a))))
}
bar_polys <- plot_data %>%
  rowwise() %>%
  mutate(poly = list(make_bar_poly(angle, inner_r, bar_len, half_w * 0.98))) %>%
  unnest(poly) %>%
  ungroup()

cls_arcs <- plot_data %>%
  group_by(Class_short) %>%
  summarise(start = min(angle) - half_w * 0.9, end = max(angle) + half_w * 0.9,
            mid = mean(angle), n = n(), .groups = "drop")

ring_r <- inner_r - 0.25
cls_ring <- do.call(rbind, lapply(seq_len(nrow(cls_arcs)), function(i) {
  a <- seq(cls_arcs$start[i], cls_arcs$end[i], length.out = 60)
  data.frame(x = ring_r * cos(a), y = ring_r * sin(a),
             Class_short = cls_arcs$Class_short[i])
}))

cls_arcs <- cls_arcs %>%
  mutate(lab_r     = ifelse(n >= 3, inner_r - 1.0, inner_r - 2.0),
         x         = lab_r * cos(mid),
         y         = lab_r * sin(mid),
         anchor_x  = (inner_r - 0.45) * cos(mid),
         anchor_y  = (inner_r - 0.45) * sin(mid),
         nudge_x   = x - anchor_x,
         nudge_y   = y - anchor_y,
         txt_angle = ifelse(cos(mid) >= 0, mid * 180 / pi, mid * 180 / pi + 180) - 90)
cls_wide   <- filter(cls_arcs, n >= 3)
cls_narrow <- filter(cls_arcs, n < 3)

lab_shift_x <- c(GlcCer = -0.45)
cls_narrow <- cls_narrow %>%
  mutate(nudge_x = nudge_x + ifelse(as.character(Class_short) %in% names(lab_shift_x),
                                    lab_shift_x[as.character(Class_short)], 0))

ring_vals <- seq(2, max_logP, by = 2)
theta     <- seq(0, 2 * pi, length.out = 361)
ref_rings <- do.call(rbind, lapply(ring_vals, function(v)
  data.frame(x = to_r(v) * cos(theta), y = to_r(v) * sin(theta), v = v)))

tick_vals <- c(0, ring_vals)
scale_df  <- data.frame(v = tick_vals, r = to_r(tick_vals))
fdr_r     <- to_r(fdr_logP)

lipid_txt_pt <- 6
class_txt_pt <- 7
scale_txt_pt <- 6.5
legend_pt    <- 7
tag_pt       <- 8

col_neg <- "#436891"
col_pos <- "#C0392B"
if (file.exists("UPFsubcat_lipids42_adj.csv")) {
  hm <- read.csv("UPFsubcat_lipids42_adj.csv", check.names = FALSE)
  beta_lim <- ceiling(max(abs(hm$beta), na.rm = TRUE) * 50) / 50
} else {
  beta_lim <- ceiling(max(abs(plot_data$beta), na.rm = TRUE) * 50) / 50
}
beta_breaks <- seq(-beta_lim, beta_lim, length.out = 5)

p_circ <- ggplot() +

  geom_path(data = ref_rings, aes(x, y, group = v),
            colour = "grey80", linewidth = 0.25, linetype = "dashed") +
  annotate("path", x = inner_r * cos(theta), y = inner_r * sin(theta),
           colour = "grey60", linewidth = 0.3) +

  annotate("rect", xmin = -inner_r * sin(gap_radians / 2) * 0.95,
           xmax = inner_r * sin(gap_radians / 2) * 0.95,
           ymin = inner_r - 0.1, ymax = inner_r + max_bar + 0.6,
           fill = "white", colour = NA) +

  geom_polygon(data = bar_polys, aes(x, y, group = lipids_name, fill = beta),
               colour = NA) +
  scale_fill_gradient2(
    low = col_neg, mid = "white", high = col_pos, midpoint = 0,
    limits = c(-beta_lim, beta_lim), breaks = beta_breaks,
    labels = gsub("-", "\u2212", sprintf("%.2f", beta_breaks)),
    oob = scales::squish, name = "Standardized \u03B2",
    guide = guide_colourbar(order = 1, direction = "horizontal",
                            title.position = "top", title.hjust = 0.5,
                            barwidth = unit(3.4, "cm"), barheight = unit(0.22, "cm"),
                            frame.colour = NA, ticks.colour = "white")) +

  geom_path(data = cls_ring, aes(x, y, group = Class_short, colour = Class_short),
            linewidth = 1.6, lineend = "butt") +
  scale_colour_manual(values = cls_colors[cls_present], breaks = cls_present,
                      labels = cls_full[cls_present], name = "Lipid class",
                      guide = guide_legend(order = 2, ncol = 2, title.position = "top",
                                           override.aes = list(linewidth = 2.5))) +

  geom_text(data = cls_wide, aes(x, y, label = Class_short, angle = txt_angle),
            colour = cls_colors[as.character(cls_wide$Class_short)],
            size = class_txt_pt / .pt, family = fig_font) +
  geom_text_repel(data = cls_narrow, aes(anchor_x, anchor_y, label = Class_short),
                  colour = cls_colors[as.character(cls_narrow$Class_short)],
                  nudge_x = cls_narrow$nudge_x, nudge_y = cls_narrow$nudge_y,
                  size = class_txt_pt / .pt, family = fig_font,
                  segment.colour = "grey50", segment.size = 0.2,
                  min.segment.length = 0, box.padding = 0.15, point.padding = 0,
                  force = 1, direction = "both", max.overlaps = Inf, seed = 41) +

  geom_text(data = plot_data,
            aes(x = (bar_len + 0.3) * cos(angle), y = (bar_len + 0.3) * sin(angle),
                label = lipids_name, angle = label_angle, hjust = hjust),
            size = lipid_txt_pt / .pt, colour = "grey10", family = fig_font) +

  annotate("segment", x = 0, xend = 0, y = inner_r, yend = to_r(max_logP),
           colour = "grey30", linewidth = 0.4) +
  geom_segment(data = scale_df, aes(x = -0.15, xend = 0.15, y = r, yend = r),
               colour = "grey30", linewidth = 0.3) +
  geom_text(data = scale_df, aes(x = -0.3, y = r, label = v), hjust = 1,
            size = scale_txt_pt / .pt, colour = "grey30", family = fig_font) +
  annotate("text", x = 0, y = to_r(max_logP) + 0.7, label = "'\u2212'*log[10]*'(P)'",
           parse = TRUE, size = scale_txt_pt / .pt, colour = "grey30", family = fig_font) +

  annotate("segment", x = -0.35, xend = 0.35, y = fdr_r, yend = fdr_r,
           colour = "black", linewidth = 0.5, linetype = "solid") +
  annotate("text", x = 0.45, y = fdr_r, label = "FDR\n0.05", hjust = 0,
           lineheight = 0.85, size = scale_txt_pt / .pt, colour = "black",
           family = fig_font) +

  geom_path(data = data.frame(
              x = fdr_r * cos(seq(pi / 2 - gap_radians / 2,
                                  pi / 2 + gap_radians / 2 - 2 * pi, length.out = 300)),
              y = fdr_r * sin(seq(pi / 2 - gap_radians / 2,
                                  pi / 2 + gap_radians / 2 - 2 * pi, length.out = 300))),
            aes(x = x, y = y),
            colour = "black", linewidth = 0.25, linetype = "dashed") +

  geom_blank(data = label_ends, aes(x, y)) +

  coord_fixed(clip = "off") +
  theme_void(base_size = 7, base_family = fig_font) +
  theme(
    legend.position   = "bottom",
    legend.box        = "vertical",
    legend.justification = "center",
    legend.box.just   = "center",
    legend.box.spacing = unit(4, "pt"),
    legend.title      = element_text(size = legend_pt, colour = "black"),
    legend.text       = element_text(size = legend_pt, colour = "black"),
    legend.key.height = unit(0.25, "cm"),
    legend.key.width  = unit(0.4, "cm"),
    legend.spacing.y  = unit(2, "pt"),
    legend.margin     = margin(0, 0, 0, 0),
    plot.tag          = element_text(size = tag_pt, face = "bold", family = fig_font),
    plot.background   = element_rect(fill = "white", colour = NA),
    plot.margin       = margin(t = 4, r = 4, b = 4, l = 4)
  )

fig_width  <- 95 / 25.4
fig_height <- 135 / 25.4

setwd(tabsFigs)

ggsave(p_circ, file = "fig_circManhattan_lipids42.tiff",
       width = fig_width, height = fig_height, units = "in", dpi = 600,
       compression = "lzw", bg = "white",
       device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_tiff else "tiff")

ggsave(p_circ, file = "fig_circManhattan_lipids42.pdf",
       width = fig_width, height = fig_height, units = "in", device = cairo_pdf)

saveRDS(list(plot = p_circ, width = fig_width, height = fig_height),
        "fig_circManhattan_lipids42.rds")

