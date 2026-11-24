# Figure 2: merge panels a, b and c.
# Input: panel .rds files from 2.1.3 and 2.1.4. Output: Figure 2 (result/tabsFigs).

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
library(tidyverse)
library(ggplot2)
library(ggtext)
library(patchwork)

fig_font <- "Arial"

result    <- file.path(proj_dir, "result")
tabsFigs  <- file.path(proj_dir, "result", "tabsFigs")

setwd(tabsFigs)
panel_a <- readRDS("fig1_panel_a_withChange.rds")
panel_b <- readRDS("fig1_panel_b_totalUPF.rds")
panel_c <- readRDS("fig1_panel_c_subclass.rds")

fig_width  <- 183 / 25.4
fig_height <- 130 / 25.4

gap_in <- 0.08

no_margin <- theme(plot.margin = margin(0, 0, 0, 0))

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
wrap_b <- add_letter(wrap_elements(full = panel_b$plot), "b")
wrap_c <- add_letter(wrap_elements(full = panel_c$plot), "c")

top_row <- (wrap_a | wrap_b) +
  plot_layout(widths = c(panel_a$width, panel_b$width))

fig1_abc <- (top_row / plot_spacer() / wrap_c) +
  plot_layout(heights = c(panel_a$height, gap_in, panel_c$height)) &
  no_margin

setwd(tabsFigs)
ggsave(
  fig1_abc,
  file        = "fig1_combined_abc.tiff",
  bg          = "white",
  width       = fig_width,
  height      = fig_height,
  units       = "in",
  dpi         = 600,
  compression = "lzw",
  device      = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_tiff else "tiff"
)

ggsave(
  fig1_abc,
  file   = "fig1_combined_abc.pdf",
  bg     = "white",
  width  = fig_width,
  height = fig_height,
  units  = "in",
  device = cairo_pdf
)
