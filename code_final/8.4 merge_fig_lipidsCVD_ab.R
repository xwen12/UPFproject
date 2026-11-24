# Figure 4: merge panels a and b.
# Input: panel .rds files from 8.1 and 8.5. Output: Figure 4 (result/tabsFigs).

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"

library(ggplot2)
library(patchwork)

tabsFigs <- file.path(proj_dir, "result", "tabsFigs")

setwd(tabsFigs)
panel_a <- readRDS("fig_lipidsCVD_panel_a_forest.rds")
panel_b_file <- "fig_lipidsCVD_panel_b_mediationPkg_survreg_main.rds"
panel_b <- readRDS(panel_b_file)

fig_w_mm <- 183
a_h_mm   <- panel_a$height * 25.4
b_w_mm   <- panel_b$width  * 25.4
b_h_mm   <- panel_b$height * 25.4
gap_mm   <- 3
fig_h_mm <- a_h_mm + gap_mm + b_h_mm

wrap_a <- wrap_elements(full = panel_a$plot)
wrap_b <- wrap_elements(full = panel_b$plot)

if (fig_w_mm - b_w_mm > 0.5) {
  bottom_row <- (wrap_b | plot_spacer()) +
    plot_layout(widths = c(b_w_mm, fig_w_mm - b_w_mm))
} else {
  bottom_row <- wrap_b
}

fig_ab <- (wrap_a / plot_spacer() / bottom_row) +
  plot_layout(heights = c(a_h_mm, gap_mm, b_h_mm)) &
  theme(plot.margin = margin(0, 0, 0, 0))

fig_w <- fig_w_mm / 25.4
fig_h <- fig_h_mm / 25.4

ggsave(fig_ab, file = "figY_lipidsCVD_ab.tiff",
       width = fig_w, height = fig_h, units = "in", dpi = 600,
       compression = "lzw", bg = "white",
       device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_tiff else "tiff")

ggsave(fig_ab, file = "figY_lipidsCVD_ab.pdf",
       width = fig_w, height = fig_h, units = "in", device = cairo_pdf)

