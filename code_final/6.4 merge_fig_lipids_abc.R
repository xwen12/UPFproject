# Figure 3: merge panels a, b and c.
# Input: panel .rds files from 6.2.1, 6.1.1 and 6.3.2. Output: Figure 3 (result/tabsFigs).

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"

library(ggplot2)
library(patchwork)

fig_font <- "Arial"

tabsFigs <- file.path(proj_dir, "result", "tabsFigs")

setwd(tabsFigs)
panel_a <- readRDS("fig_circManhattan_lipids42.rds")
panel_b <- readRDS("fig_classBeta_lipids.rds")
panel_c <- readRDS("fig_heatmap_lipids42_UPFsubcat.rds")

for (nm in c("panel_a", "panel_b", "panel_c")) {
  p <- get(nm)
}

plot_b <- panel_b$plot

plot_c <- panel_c$plot

panel_letter <- function(lab) {
  grid::textGrob(lab, x = unit(1.5, "mm"), y = unit(1, "npc") - unit(1.5, "mm"),
                 hjust = 0, vjust = 1,
                 gp = grid::gpar(fontsize = 8, fontface = "bold", fontfamily = fig_font))
}
add_letter <- function(p, lab) {
  p + inset_element(wrap_elements(full = panel_letter(lab)) +
                      theme(plot.background = element_blank()),
                    left = 0, bottom = 0, right = 1, top = 1,
                    align_to = "full", on_top = TRUE, ignore_tag = TRUE)
}

wrap_a <- add_letter(wrap_elements(full = panel_a$plot), "a")
wrap_b <- add_letter(wrap_elements(full = plot_b), "b")
wrap_c <- add_letter(wrap_elements(full = plot_c), "c")

fig_w_mm <- 183
fig_h_mm <- 235

top_h_mm <- 152
a_w_mm   <- 95

top_row <- (wrap_a | wrap_b) +
  plot_layout(widths = c(a_w_mm, fig_w_mm - a_w_mm))

gap_mm  <- 5
fig_abc <- (top_row / plot_spacer() / wrap_c) +
  plot_layout(heights = c(top_h_mm, gap_mm, fig_h_mm - top_h_mm - gap_mm)) &
  theme(plot.margin = margin(1, 1, 1, 1))

fig_w <- fig_w_mm / 25.4
fig_h <- fig_h_mm / 25.4

ggsave(fig_abc, file = "figX_lipids_abc.tiff",
       width = fig_w, height = fig_h, units = "in", dpi = 600,
       compression = "lzw", bg = "white",
       device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_tiff else "tiff")

ggsave(fig_abc, file = "figX_lipids_abc.pdf",
       width = fig_w, height = fig_h, units = "in", device = cairo_pdf)

