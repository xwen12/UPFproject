# Figure 3c: heatmap of betas for total UPF and 8 UPF subclasses with the 42 UPF-related lipids.
# Input: 6.3.1 results, lipid class annotation. Output: panel c (result/tabsFigs).

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
raw_dir  <- "path/to/SHS_raw_data"

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)

fig_font <- "Arial"

result   <- file.path(proj_dir, "result")
tabsFigs <- file.path(proj_dir, "result", "tabsFigs")
lipid_info_path <- file.path(raw_dir, "lipidomics data", "lipid_4_6_compare.csv")

setwd(result)
res <- read.csv("UPFsubcat_lipids42_adj.csv", header = TRUE, stringsAsFactors = FALSE,
                check.names = FALSE)

nova_short <- c(
  "nova4_kcal_P"   = "Total UPF",
  "nova4_7_kcal_P" = "Ready-to-eat mixed dishes",
  "nova4_3_kcal_P" = "Sweet snacks & desserts",
  "nova4_5_kcal_P" = "Sweetened beverages",
  "nova4_2_kcal_P" = "Sauces & spreads",
  "nova4_6_kcal_P" = "Ready-to-eat meat products",
  "nova4_1_kcal_P" = "Breads & breakfast foods",
  "nova4_4_kcal_P" = "Savory snacks",
  "nova4_8_kcal_P" = "Dairy desserts"
)

nova_colors <- c(
  "Total UPF"                  = "#1A1A1A",
  "Ready-to-eat mixed dishes"  = "#A7799A",
  "Sweet snacks & desserts"    = "#E0BE4F",
  "Sweetened beverages"        = "#74AAA5",
  "Sauces & spreads"           = "#DD5454",
  "Ready-to-eat meat products" = "#5FA64C",
  "Breads & breakfast foods"   = "#F29533",
  "Savory snacks"              = "#436891",
  "Dairy desserts"             = "#AAAAAA"
)

res <- res %>%
  mutate(nova_short = factor(unname(nova_short[nova]), levels = unname(nova_short)))

lipid_info <- read.csv(lipid_info_path, header = TRUE, stringsAsFactors = FALSE,
                       check.names = FALSE)

lipid_meta <- data.frame(lipids_name = unique(res$lipids_name), stringsAsFactors = FALSE) %>%
  left_join(lipid_info[, c("lipids_name", "Class_short")], by = "lipids_name")

lipid_meta[lipid_meta$Class_short == "N/A",                       "Class_short"] <- "Unknown"
lipid_meta[lipid_meta$Class_short %in% c("FAHFA", "scFA", "USFA"), "Class_short"] <- "FA"
lipid_meta[lipid_meta$Class_short %in% c("Glc"),                  "Class_short"] <- "GlcCer"

cls_order <- c("AC", "CE", "CER", "DAG", "FA", "GlcCer", "PC", "PE", "PG", "PI",
               "PS", "SM", "TAG", "Unknown")
cls_order <- cls_order[cls_order %in% lipid_meta$Class_short]

lipid_meta <- lipid_meta %>%
  mutate(Class_short = factor(Class_short, levels = cls_order)) %>%
  arrange(Class_short, lipids_name)

gap <- 0.6
lipid_meta <- lipid_meta %>%
  mutate(cls_id = as.integer(Class_short),
         x      = row_number() + (cls_id - 1) * gap)

n_exp   <- nlevels(res$nova_short)
exp_pos <- data.frame(nova_short = factor(levels(res$nova_short), levels = levels(res$nova_short)),
                      y = rev(seq_len(n_exp)))

plot_df <- res %>%
  left_join(lipid_meta %>% select(lipids_name, Class_short, x), by = "lipids_name") %>%
  left_join(exp_pos, by = "nova_short")

lipid_txt_pt <- 6
axis_txt_pt  <- 7
class_txt_pt <- 6.5
sig_txt_pt   <- 6
legend_pt    <- 7
tag_pt       <- 8

col_neg <- "#436891"
col_pos <- "#C0392B"

beta_lim    <- ceiling(max(abs(plot_df$beta), na.rm = TRUE) * 50) / 50
beta_breaks <- seq(-beta_lim, beta_lim, length.out = 5)

cls_colors <- c(
  AC = "#3F7F33", CE = "#2F4F73", CER = "#C94C7C", DAG = "#4C78A8",
  FA = "#B8860B", GlcCer = "#4D4D4D", PC = "#2C6DA4", PE = "#6A5AAE",
  PG = "#C75D1A", PI = "#AD3F7A", PS = "#8E6BBF", SM = "#1F7A6E",
  TAG = "#A0522D", Unknown = "#7A7A7A")

band_y <- n_exp + 0.95
cls_band <- lipid_meta %>%
  group_by(Class_short) %>%
  summarise(xmin = min(x) - 0.5, xmax = max(x) + 0.5, xmid = mean(x),
            n = n(), .groups = "drop") %>%
  mutate(col = cls_colors[as.character(Class_short)])
cls_in  <- filter(cls_band, n >= 3)
cls_out <- filter(cls_band, n <  3)

key_x  <- min(lipid_meta$x) - 1.3
key_df <- exp_pos %>% mutate(col = nova_colors[as.character(nova_short)])

heatmap_lipids <- ggplot(plot_df, aes(x = x, y = y)) +

  geom_tile(aes(fill = beta), colour = "white", linewidth = 0.25) +

  geom_text(aes(label = sig), size = sig_txt_pt / .pt, vjust = 0.75,
            colour = "black", family = fig_font) +

  annotate("segment", x = min(lipid_meta$x) - 0.5, xend = max(lipid_meta$x) + 0.5,
           y = n_exp - 0.5, yend = n_exp - 0.5, colour = "black", linewidth = 0.4) +

  annotate("rect", xmin = cls_band$xmin, xmax = cls_band$xmax,
           ymin = band_y - 0.35, ymax = band_y + 0.35,
           fill = cls_band$col, colour = NA) +
  annotate("text", x = cls_in$xmid, y = band_y, label = as.character(cls_in$Class_short),
           colour = "white", size = class_txt_pt / .pt, family = fig_font) +
  annotate("text", x = cls_out$xmid, y = band_y + 0.5,
           label = as.character(cls_out$Class_short), colour = cls_out$col,
           angle = 45, hjust = 0, vjust = 0, size = class_txt_pt / .pt, family = fig_font) +

  annotate("tile", x = key_x, y = key_df$y, width = 0.9, height = 0.9,
           fill = key_df$col) +

  scale_fill_gradient2(
    low = col_neg, mid = "white", high = col_pos, midpoint = 0,
    limits = c(-beta_lim, beta_lim), breaks = beta_breaks,
    labels = gsub("-", "\u2212", sprintf("%.2f", beta_breaks)),
    oob = scales::squish, na.value = "grey85",
    name = "Standardized\n\u03B2",
    guide = guide_colourbar(direction = "vertical", title.position = "top",
                            barwidth = unit(0.22, "cm"), barheight = unit(2.0, "cm"),
                            frame.colour = NA, ticks.colour = "white")) +

  scale_x_continuous(breaks = lipid_meta$x, labels = lipid_meta$lipids_name,
                     expand = c(0, 0.3)) +
  scale_y_continuous(breaks = exp_pos$y, labels = as.character(exp_pos$nova_short),
                     expand = c(0, 0.1)) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +

  theme_classic(base_size = 7, base_family = fig_font) +
  theme(
    axis.line          = element_blank(),
    axis.ticks         = element_blank(),
    axis.text.x        = element_text(size = lipid_txt_pt, colour = "black",
                                      angle = 90, hjust = 1, vjust = 0.5,
                                      margin = margin(t = 1)),
    axis.text.y        = element_text(size = axis_txt_pt, colour = "black",
                                      margin = margin(r = 1)),
    legend.position    = "right",
    legend.title       = element_text(size = legend_pt, colour = "black"),
    legend.text        = element_text(size = legend_pt, colour = "black"),
    legend.box.spacing = unit(4, "pt"),
    plot.tag           = element_text(size = tag_pt, face = "bold", family = fig_font),
    plot.background    = element_rect(fill = "white", colour = NA),
    panel.background   = element_rect(fill = "white", colour = NA),
    plot.margin        = margin(t = 4, r = 4, b = 4, l = 4)
  )

g_hm    <- ggplotGrob(heatmap_lipids)
leg_idx <- grep("guide-box", g_hm$layout$name)
leg_idx <- leg_idx[!sapply(g_hm$grobs[leg_idx], inherits, "zeroGrob")][1]
key_grob  <- g_hm$grobs[[leg_idx]]
note_grob <- grid::textGrob("** FDR < 0.05\n* P < 0.05", x = unit(1, "mm"), y = 1,
                            hjust = 0, vjust = 1,
                            gp = grid::gpar(fontsize = legend_pt, fontfamily = fig_font,
                                            lineheight = 0.95))
key_h_mm  <- grid::convertHeight(grid::grobHeight(key_grob), "mm", valueOnly = TRUE)
key_top_mm <- 12
key_block <- gtable::gtable_matrix(
  "key_block",
  matrix(list(grid::nullGrob(), key_grob, grid::nullGrob(), note_grob, grid::nullGrob()),
         ncol = 1),
  widths  = unit(1, "null"),
  heights = grid::unit.c(unit(c(key_top_mm, key_h_mm, 3, 7), "mm"),
                         unit(1, "null")),
  clip    = "off")

fig_heatmap <- (heatmap_lipids + theme(legend.position = "none")) +
  wrap_elements(full = key_block) +
  plot_layout(widths = c(1, 0.10))

fig_width  <- 183 / 25.4
fig_height <- 80 / 25.4

setwd(tabsFigs)

ggsave(fig_heatmap, file = "fig_heatmap_lipids42_UPFsubcat.tiff",
       width = fig_width, height = fig_height, units = "in", dpi = 600,
       compression = "lzw", bg = "white",
       device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_tiff else "tiff")

ggsave(fig_heatmap, file = "fig_heatmap_lipids42_UPFsubcat.pdf",
       width = fig_width, height = fig_height, units = "in", device = cairo_pdf)

saveRDS(list(plot = fig_heatmap, width = fig_width, height = fig_height),
        "fig_heatmap_lipids42_UPFsubcat.rds")

