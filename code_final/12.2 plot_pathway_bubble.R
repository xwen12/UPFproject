# Figure 5: bubble plot of KEGG pathway results (MetaboAnalyst 6.0 and LIPEA).
# Input: exported pathway results (result/enrich). Output: result/tabsFigs/pathway_bubble_MA_LIPEA.*

rm(list = ls())

proj_dir <- "path/to/proj16_UPF_rework"
library(dplyr)
library(stringr)
library(ggplot2)
library(readxl)
library(patchwork)
library(scales)

tabsFigs      <- file.path(proj_dir, "result", "tabsFigs")
enrich_output <- file.path(proj_dir, "result", "enrich")

if (!dir.exists(tabsFigs)) dir.create(tabsFigs, recursive = TRUE)

setwd(enrich_output)

ma_raw <- read.csv('0metaboanalyst_pathway_results.csv', header = TRUE,
                   check.names = FALSE, stringsAsFactors = FALSE)
names(ma_raw)[1] <- "Pathway"

lipea_raw <- read_excel('0lipea_results.xls', col_names = FALSE, .name_repair = "minimal")
hdr_row   <- which(lipea_raw[[1]] == "Pathway name")[1]
names(lipea_raw) <- as.character(unlist(lipea_raw[hdr_row, ]))
lipea_raw <- lipea_raw[-(1:hdr_row), ] %>%
  filter(!is.na(`Pathway name`))

ma <- ma_raw %>%
  transmute(tool    = "MetaboAnalyst",
            Pathway = Pathway,
            hits    = as.numeric(Hits),
            total   = as.numeric(Total),
            p       = as.numeric(`Raw p`),
            FDR     = as.numeric(FDR),
            impact  = as.numeric(Impact))

lipea <- lipea_raw %>%
  transmute(tool     = "LIPEA",
            Pathway  = `Pathway name`,
            hits     = as.numeric(`Converted lipids (number)`),
            total    = as.numeric(`Pathway lipids`),
            p        = as.numeric(`p-value`),
            FDR      = as.numeric(`Benjamini correction`),
            impact   = NA_real_,
            kegg_ids = `Converted lipids (list)`)

ma    <- ma    %>% mutate(ratio = hits / total)
lipea <- lipea %>% mutate(ratio = hits / total)

ma_plot <- ma

lipea_drop <- c("Pathogenic Escherichia coli infection",
                "Choline metabolism in cancer",
                "Autophagy - other",
                "Leishmaniasis",
                "Tuberculosis",
                "Kaposi's sarcoma-associated herpesvirus infection")
lipea_plot <- lipea %>%
  filter(p < 0.05, !Pathway %in% lipea_drop)

short_name <- function(x) {
  x %>%
    str_replace("^Glycosylphosphatidylinositol \\(GPI\\)-anchor biosynthesis$",
                "GPI-anchor biosynthesis") %>%
    str_replace("^alpha-Linolenic", "\u03b1-Linolenic")
}

prep_plot <- function(d) {
  d %>%
    mutate(Pathway_lab = short_name(Pathway),
           logp        = -log10(p),
           logFDR      = -log10(FDR),
           fdr_grp     = factor(ifelse(FDR < 0.05, "sig", "ns"),
                                levels = c("sig", "ns"))) %>%
    arrange(p) %>%
    mutate(Pathway_lab = factor(Pathway_lab, levels = rev(Pathway_lab)))
}
ma_plot    <- prep_plot(ma_plot)
lipea_plot <- prep_plot(lipea_plot)

plot_data <- bind_rows(ma_plot, lipea_plot) %>%
  select(tool, Pathway, hits, total, ratio, impact, p, FDR, kegg_ids)
write.table(plot_data, 'pathway_bubble_plotdata.csv', sep = ',', row.names = F, col.names = T)

col_low  <- "#436891"
col_mid  <- "#B5B5B5"
col_hi0  <- "#E8A397"
col_hi1  <- "#D9604F"
col_high <- "#C0392B"

x_max   <- max(c(ma_plot$logp, lipea_plot$logp))
col_lim <- c(0, ceiling(max(c(ma_plot$logFDR, lipea_plot$logFDR))))
col_val <- rescale(c(0, -log10(0.05), 1.5, 2.5, col_lim[2]), from = col_lim)

x_lim <- c(0, x_max * 1.08)

fdr_shapes <- c("sig" = 16, "ns" = 1)
fdr_labels <- c(expression("FDR" < 0.05), expression("FDR" >= 0.05))

theme_nat <- theme_classic(base_size = 9, base_family = "sans") +
  theme(
    panel.background   = element_rect(fill = "white", colour = NA),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.line          = element_line(colour = "grey80", linewidth = 0.5),
    axis.ticks         = element_blank(),
    axis.text          = element_text(size = 8, colour = "black"),
    axis.title.x       = element_text(size = 8.5, colour = "black", margin = margin(t = 4)),
    plot.title         = element_text(size = 9, face = "bold", hjust = 0),
    plot.tag           = element_text(size = 11, face = "bold", family = "sans"),
    legend.position    = "right",
    legend.title       = element_text(size = 7.5),
    legend.text        = element_text(size = 7),
    legend.key.size    = unit(0.35, "cm"),
    legend.spacing.y   = unit(0.1, "cm"),
    plot.margin        = margin(t = 4, r = 6, b = 4, l = 4)
  )

plot_bubble <- function(d, title_txt, size_var, size_name, size_lim, size_brk,
                        size_lab = waiver()) {
  ggplot(d, aes(x = logp, y = Pathway_lab)) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed",
               colour = "grey55", linewidth = 0.35) +
    geom_point(aes(size = .data[[size_var]], colour = logFDR, shape = fdr_grp),
               stroke = 0.9) +
    scale_x_continuous(limits = x_lim, expand = expansion(mult = c(0, 0.02)),
                       breaks = pretty_breaks(n = 5)) +
    scale_size(name = size_name, limits = size_lim, breaks = size_brk,
               labels = size_lab, range = c(1.3, 6)) +
    scale_colour_gradientn(name = expression(-log[10]*"(FDR)"),
                           colours = c(col_low, col_mid, col_hi0, col_hi1, col_high),
                           values = col_val, limits = col_lim,
                           breaks = pretty_breaks(n = 4)) +
    scale_shape_manual(name = NULL, values = fdr_shapes, labels = fdr_labels,
                       drop = FALSE) +
    labs(x = expression(-log[10]*"("*italic(P)*")"), y = NULL, title = title_txt) +
    guides(shape  = guide_legend(order = 1,
                                 override.aes = list(size = 2.5, colour = "grey40")),
           colour = guide_colourbar(order = 2, barwidth = unit(0.3, "cm"),
                                    barheight = unit(1.8, "cm")),
           size   = guide_legend(order = 3,
                                 override.aes = list(shape = 16, colour = "grey60"))) +
    theme_nat
}

imp_max <- max(ma_plot$impact)
p_ma <- plot_bubble(ma_plot,
                    title_txt = "MetaboAnalyst pathway enrichment",
                    size_var  = "impact",
                    size_name = "Pathway impact",
                    size_lim  = c(0, imp_max),
                    size_brk  = pretty(c(0, imp_max), n = 3))

p_lipea <- plot_bubble(lipea_plot,
                       title_txt = "LIPEA pathway enrichment",
                       size_var  = "ratio",
                       size_name = "Enrichment ratio",
                       size_lim  = c(0, max(lipea_plot$ratio)),
                       size_brk  = setdiff(pretty(c(0, max(lipea_plot$ratio)), n = 3), 0),
                       size_lab  = percent_format(accuracy = 1))

n_a <- nrow(ma_plot)
n_b <- nrow(lipea_plot)
fig_pathway <- (p_ma / p_lipea) +
  plot_layout(heights = c(n_a + 1, n_b + 1)) +
  plot_annotation(tag_levels = "a")

fig_width  <- 5.8
fig_height <- 1.6 + 0.2 * (n_a + n_b + 2)

setwd(tabsFigs)
ggsave(fig_pathway, file = "pathway_bubble_MA_LIPEA.tiff", bg = "white",
       width = fig_width, height = fig_height, units = "in",
       dpi = 350, compression = "lzw")
ggsave(fig_pathway, file = "pathway_bubble_MA_LIPEA.pdf", device = cairo_pdf,
       width = fig_width, height = fig_height, units = "in")
