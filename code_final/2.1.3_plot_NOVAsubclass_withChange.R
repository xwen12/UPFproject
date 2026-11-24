# Figure 2a: mean % energy from total UPF and 8 UPF subclasses at baseline and follow-up.
# Run order: 2.1.3 -> 2.1.4 -> 2.1.5. Output: panel a (result/tabsFigs).

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
raw_dir  <- "path/to/SHS_raw_data"
library(tidyverse)
library(ggplot2)
library(ggtext)
library(patchwork)

fig_font <- "Arial"

rawData   <- raw_dir
readyData <- file.path(proj_dir, "data")
result    <- file.path(proj_dir, "result")
tabsFigs  <- file.path(proj_dir, "result", "tabsFigs")

if (!dir.exists(tabsFigs)) dir.create(tabsFigs, recursive = TRUE)

setwd(readyData)
load(file = "dat_AimOne_p4.Rda")
load(file = "dat_AimOne_p5.Rda")

cols <- paste0("nova4_", 1:8, "_kcal_P")

get_means <- function(data, group_name) {
  mean_vals <- data %>%
    summarise(across(all_of(cols), ~mean(.x, na.rm = TRUE)))
  tibble(
    category   = cols,
    mean_value = as.numeric(mean_vals) * 100,
    group      = group_name
  )
}
combined_mean <- bind_rows(
  get_means(dat_AimOne_p4, "p4"),
  get_means(dat_AimOne_p5, "p5")
)

category_labels <- c(
  "nova4_1_kcal_P" = "Ultra-processed breads & breakfast foods",
  "nova4_2_kcal_P" = "Sauces, spreads & condiments",
  "nova4_3_kcal_P" = "Packaged sweet snacks & desserts",
  "nova4_4_kcal_P" = "Packaged savory snacks",
  "nova4_5_kcal_P" = "Sugar/artificially sweetened beverages",
  "nova4_6_kcal_P" = "Meat/poultry/seafood-based ready-to-eat products",
  "nova4_7_kcal_P" = "Ready-to-eat/heat mixed dishes",
  "nova4_8_kcal_P" = "Dairy-based desserts"
)
combined_mean <- combined_mean %>%
  mutate(category = recode(category, !!!category_labels))

stack_order <- c(
  "Dairy-based desserts",
  "Packaged savory snacks",
  "Ultra-processed breads & breakfast foods",
  "Meat/poultry/seafood-based ready-to-eat products",
  "Sauces, spreads & condiments",
  "Sugar/artificially sweetened beverages",
  "Packaged sweet snacks & desserts",
  "Ready-to-eat/heat mixed dishes"
)

combined_mean <- combined_mean %>%
  mutate(
    category = factor(category, levels = rev(stack_order)),
    group    = factor(group,
                      levels = c("p4", "p5"),
                      labels = c("Baseline\n(2001–2003)", "Follow-up\n(2006–2009)"))
  )

label_df <- combined_mean %>%
  arrange(group, desc(as.integer(category))) %>%
  group_by(group) %>%
  mutate(
    cum_top = cumsum(mean_value),
    cum_bot = cum_top - mean_value,
    mid_y   = (cum_top + cum_bot) / 2
  ) %>%
  ungroup()

large_labels <- label_df %>% filter(mean_value >= 3)
small_labels <- label_df %>% filter(mean_value < 3)

bar_half  <- 0.55 / 2
small_gap <- 4.0
small_dx  <- 0.12
small_y0  <- 1.2

spread_y <- function(y, gap, y0) {
  out <- y
  out[1] <- max(y[1], y0)
  if (length(y) > 1) {
    for (i in 2:length(y)) out[i] <- max(y[i], out[i - 1] + gap)
  }
  out
}

small_labels <- small_labels %>%
  arrange(group, mid_y) %>%
  group_by(group) %>%
  mutate(label_y = spread_y(mid_y, small_gap, small_y0)) %>%
  ungroup() %>%
  mutate(
    bar_x      = as.numeric(group),
    seg_x_from = ifelse(bar_x == 1,
                        bar_x - bar_half,
                        bar_x + bar_half),
    label_x    = ifelse(bar_x == 1,
                        bar_x - bar_half - small_dx,
                        bar_x + bar_half + small_dx),
    line_x_end = ifelse(bar_x == 1,
                        label_x + 0.02,
                        label_x - 0.02),
    hjust_val  = ifelse(bar_x == 1, 1, 0)
  )

totals <- combined_mean %>%
  group_by(group) %>%
  summarise(total = sum(mean_value), .groups = "drop")

total_delta <- round(totals$total[2], 1) - round(totals$total[1], 1)
total_delta_lab <- paste0(ifelse(total_delta >= 0, "+", ""),
                          formatC(total_delta, digits = 1, format = "f"), "%")

delta_df <- combined_mean %>%
  select(category, group, mean_value) %>%
  pivot_wider(names_from = group, values_from = mean_value) %>%
  rename(ph4 = `Baseline\n(2001–2003)`,
         ph5 = `Follow-up\n(2006–2009)`) %>%
  mutate(
    delta     = round(ph5, 1) - round(ph4, 1),
    direction = ifelse(delta >= 0, "Increase", "Decrease"),
    delta_lab = ifelse(delta >= 0,
                       paste0("+", formatC(delta, digits = 1, format = "f")),
                       formatC(delta, digits = 1, format = "f"))
  )

nova4_colors <- c(
  "Packaged sweet snacks & desserts"                 = "#E0BE4F",
  "Sugar/artificially sweetened beverages"           = "#74AAA5",
  "Dairy-based desserts"                             = "#AAAAAA",
  "Ultra-processed breads & breakfast foods"         = "#F29533",
  "Sauces, spreads & condiments"                     = "#DD5454",
  "Packaged savory snacks"                           = "#436891",
  "Meat/poultry/seafood-based ready-to-eat products" = "#5FA64C",
  "Ready-to-eat/heat mixed dishes"                   = "#A7799A"
)

bar_num_pt    <- 7
total_pt      <- 7
axis_txt_pt   <- 7
axis_title_pt <- 7
tag_pt        <- 8

direction_colors <- c("Increase" = "#C0392B", "Decrease" = "#436891")

y_top    <- ceiling(max(totals$total) + 6)
y_breaks <- seq(0, floor(max(totals$total) / 10) * 10, 10)

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

subclass_barplot <- ggplot(
  combined_mean,
  aes(x = group, y = mean_value, fill = category)
) +

  geom_col(position = "stack", width = 0.55,
           colour = "white", linewidth = 0) +

  geom_text(
    data        = large_labels,
    aes(x       = group,
        y       = mid_y,
        label   = formatC(mean_value, digits = 1, format = "f")),
    inherit.aes = FALSE,
    size        = bar_num_pt / .pt,
    color       = "white",
    fontface    = "bold",
    family      = fig_font
  ) +

  geom_segment(
    data        = small_labels,
    aes(x    = seg_x_from, xend = line_x_end,
        y    = mid_y,      yend = label_y),
    inherit.aes = FALSE,
    color       = "grey40",
    linewidth   = 0.25
  ) +

  geom_text(
    data        = small_labels,
    aes(x     = label_x,
        y     = label_y,
        label = formatC(mean_value, digits = 1, format = "f"),
        hjust = hjust_val),
    inherit.aes = FALSE,
    size        = bar_num_pt / .pt,
    color       = "grey20",
    fontface    = "bold",
    family      = fig_font
  ) +

  geom_text(
    data        = totals,
    aes(x = group, y = total,
        label = paste0("Total: ", formatC(total, digits = 1, format = "f"), "%")),
    inherit.aes = FALSE,
    vjust       = -0.6,
    size        = total_pt / .pt,
    fontface    = "bold",
    family      = fig_font,
    color       = "black"
  ) +

  geom_segment(
    data        = totals,
    aes(x    = as.numeric(group) - 0.35,
        xend = as.numeric(group) + 0.35,
        y    = total, yend = total),
    inherit.aes = FALSE,
    linetype    = "dashed",
    color       = "black",
    linewidth   = 0.3
  ) +

  scale_fill_manual(values = nova4_colors, guide = "none") +

  scale_x_discrete(expand = expansion(add = c(0.6, 0.6))) +

  scale_y_continuous(
    limits = c(0, y_top),
    breaks = y_breaks,
    labels = y_breaks,
    expand = expansion(mult = c(0, 0))
  ) +

  labs(
    x    = NULL,
    y    = "% Energy from UPF",
    tag  = NULL
  ) +

  theme_classic(base_size = 7, base_family = fig_font) +
  theme(
    panel.background   = element_rect(fill = "white", colour = NA),
    panel.border       = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.line          = element_line(colour = "grey80", linewidth = 0.3),
    axis.ticks.x       = element_blank(),
    axis.ticks.y       = element_line(colour = "grey80", linewidth = 0.3),
    axis.ticks.length  = unit(0.08, "cm"),
    axis.text          = element_text(size = axis_txt_pt, colour = "black"),
    axis.title.y       = element_text(size = axis_title_pt, colour = "black", margin = margin(r = 3)),
    plot.tag           = element_text(size = tag_pt, face = "bold", family = fig_font),
    plot.tag.position  = "topleft",
    plot.margin        = margin(t = 2, r = 0, b = 4, l = 2)
  )

fig_height    <- 2.45
fig_width     <- 4.35
legend_txt_pt <- 7
header_pt     <- 7
note_pt       <- 6
row_in        <- 0.19
key_x         <- 0.05
name_x        <- 0.13
gap_in        <- 0.06
right_pad_in  <- 0.06

legend_df <- delta_df %>%
  mutate(category = factor(as.character(category), levels = rev(stack_order))) %>%
  arrange(category) %>%
  mutate(
    row_id    = row_number(),
    name_txt  = recode(as.character(category), !!!short_labels),
    delta_txt = paste0("(", str_replace(delta_lab, "^-", "−"), "%)")
  )

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
  ragg::agg_tiff(tmp_tiff, width = 5, height = 5, units = "in", res = 350)
} else {
  tiff(tmp_tiff, width = 5, height = 5, units = "in", res = 350)
}
max_name_in  <- 1.03 * max(text_width_in(legend_df$name_txt, legend_txt_pt))
max_delta_in <- 1.03 * max(text_width_in(legend_df$delta_txt, legend_txt_pt, face = "bold"))
header_in    <- 1.03 * (text_width_in("UPF subclass ", header_pt, face = "bold") +
                        text_width_in("(\u0394 % energy, baseline to follow-up)", note_pt))
invisible(dev.off())
unlink(tmp_tiff)

delta_x   <- name_x + max_name_in + gap_in
legend_w  <- max(delta_x + max_delta_in, header_in) + right_pad_in

panel_h   <- fig_height - (2 + 4) / 72
n_rows    <- nrow(legend_df)
block_h   <- (n_rows - 1) * row_in + 1.2 * row_in
top_y     <- panel_h / 2 + block_h / 2
legend_df <- legend_df %>%
  mutate(y_pos = top_y - 1.2 * row_in - (row_id - 1) * row_in)

header_lab <- paste0(
  "<span style='font-size:", header_pt, "pt'><b>UPF subclass</b></span> ",
  "<span style='font-size:", note_pt, "pt; color:#4D4D4D'>",
  "(&#916; % energy, baseline to follow-up)</span>"
)

legend_panel <- ggplot(legend_df) +

  geom_point(
    aes(x = key_x, y = y_pos, fill = category),
    shape  = 22,
    size   = 2.3,
    colour = "white",
    stroke = 0
  ) +

  geom_text(
    aes(x = name_x, y = y_pos, label = name_txt),
    hjust  = 0,
    size   = legend_txt_pt / .pt,
    colour = "grey15",
    family = fig_font
  ) +

  geom_text(
    aes(x = delta_x, y = y_pos, label = delta_txt),
    colour   = "grey15",
    hjust    = 0,
    size     = legend_txt_pt / .pt,
    fontface = "bold",
    family   = fig_font
  ) +

  ggtext::geom_richtext(
    data          = data.frame(x = 0, y = top_y, lab = header_lab),
    aes(x = x, y = y, label = lab),
    hjust         = 0,
    family        = fig_font,
    fill          = NA,
    label.color   = NA,
    label.padding = unit(0, "pt")
  ) +

  scale_fill_manual(values = nova4_colors, guide = "none") +
  scale_x_continuous(limits = c(0, legend_w), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, panel_h),  expand = c(0, 0)) +
  coord_cartesian(clip = "off") +

  theme_void(base_family = fig_font) +
  theme(plot.margin = margin(t = 2, r = 0, b = 4, l = 0))

fig1_combined <- subclass_barplot + legend_panel +
  plot_layout(widths = unit(c(1, legend_w), c("null", "in")))

setwd(tabsFigs)
ggsave(
  fig1_combined,
  file        = "fig1_panel_a_withChange.tiff",
  bg          = "white",
  width       = fig_width,
  height      = fig_height,
  units       = "in",
  dpi         = 600,
  compression = "lzw",
  device      = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_tiff else "tiff"
)

saveRDS(list(plot = patchwork::patchworkGrob(fig1_combined), width = fig_width, height = fig_height),
        file = "fig1_panel_a_withChange.rds")
