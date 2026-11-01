source(here::here("Code", "pipeline_paths.R"))
library(ggplot2)
library(ggpubr)
library(scales)
library(here)

load(proteinfix_input('met_a26_PCAdata.RData'))
load(proteinfix_input('met_a27_PCAdata.RData'))
load(proteinfix_input('pro_a26_PCAdata.RData'))
load(proteinfix_input('pro_a27_PCAdata.RData'))

# -----------------------------
# Figure export parameters (Fig 2 = post, S2 = pre)
# -----------------------------
## ---- shared PCA style: keep IDENTICAL across Fig 2 / Fig 5 / S2 ----
pca_pt_size    <- 2    # in-plot point size (ggplot units)
pca_legend_pt  <- 3    # legend key point size (matches existing 4-panel)
pca_legend_txt <- 9    # legend text size, pt (matches existing 4-panel)

# SR: single uploaded file per figure, TIFF/LZW/RGB, >=300 dpi at final size.
# Export at intended placement width so pt-sized text lands at true size and
# 600 dpi is 600 dpi *there*. 180 mm = full-width (verify against SR spec).
fig_fmt  <- "tiff"
fig_dpi  <- 600
fig_w_mm <- 183      # double-column width (your spec); was 250 mm, SR would have downscaled
fig_h_mm <- 146      # holds ~1.25 aspect at 183 mm (183/1.25=146.4)

# Panel spacing: extra vertical margin per panel so the top row's x-axis title
# doesn't crowd the bottom row. Applied at assembly (make_pca_plot untouched).
# Tunable — render and adjust this one number; default plot.margin is ~5.5 pt.
panel_gap_v_pt <- 12   # top/bottom room (already good at 12)
panel_gap_h_pt <- 12   # left/right room: clears right column's y-title AND stops the
                       # far-right x tick label (e.g. "200" on the pre/S2 figure) clipping.
                       # Was 4 pt (below ggplot's 5.5 default), which caused both problems.
panel_margin <- theme(plot.margin = margin(t = panel_gap_v_pt, r = panel_gap_h_pt,
                                            b = panel_gap_v_pt, l = panel_gap_h_pt))
output_dir <- getOption("proteinfix.pca_export_dir", here::here("Figure"))   # conform to repo convention; run from project root so this resolves to <root>/Figure
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Shared legend levels/colors/shapes
# -----------------------------
legend_levels <- c(
  paste("Plate", 1:7),
  "SPQC"
)

legend_colors <- c(
  "Plate 1" = "#E69F00",
  "Plate 2" = "#56B4E9",
  "Plate 3" = "#009E73",
  "Plate 4" = "#F0E442",
  "Plate 5" = "#0072B2",
  "Plate 6" = "#D55E00",
  "Plate 7" = "#CC79A7",
  "SPQC"    = "#333333"
)

legend_shapes <- c(
  "Plate 1" = 16,
  "Plate 2" = 16,
  "Plate 3" = 16,
  "Plate 4" = 16,
  "Plate 5" = 16,
  "Plate 6" = 16,
  "Plate 7" = 16,
  "SPQC"    = 17
)

# Match the manuscript layout: shared legend, column headings, and row labels.
assemble_pca_figure <- function(legend_plot, panels) {
  cowplot::ggdraw() +
    cowplot::draw_plot(legend_plot, x = 0.06, y = 0.94, width = 0.94, height = 0.06) +
    cowplot::draw_label("Metabolomic", x = 0.295, y = 0.915, size = 14) +
    cowplot::draw_label("Proteomic", x = 0.765, y = 0.915, size = 14) +
    cowplot::draw_plot(panels, x = 0.06, y = 0, width = 0.94, height = 0.89) +
    cowplot::draw_label("A26", x = 0.028, y = 0.6675, size = 14) +
    cowplot::draw_label("A27", x = 0.028, y = 0.2225, size = 14)
}

################################ Figures for Post-Combat ################################################

# -----------------------------
# Common plotting function
# -----------------------------
make_pca_plot <- function(data_mat, legend_group, title_text) {
  pca_res <- prcomp(data_mat, center = TRUE, scale. = TRUE)
  pca_df <- as.data.frame(pca_res$x[, 1:2])
  colnames(pca_df) <- c("PC1", "PC2")

  var_exp <- 100 * (pca_res$sdev^2 / sum(pca_res$sdev^2))
  message(title_text, ": n=", nrow(data_mat), ", analytes=", ncol(data_mat),
          ", PC1=", round(var_exp[1], 1), "%, PC2=", round(var_exp[2], 1), "%")

  pca_df$LegendGroup <- factor(legend_group, levels = legend_levels)

  ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = LegendGroup, shape = LegendGroup), size = pca_pt_size) +
    scale_color_manual(values = legend_colors, drop = FALSE) +
    scale_shape_manual(values = legend_shapes, drop = FALSE) +
    scale_x_continuous(breaks = pretty_breaks(n = 4)) +
    scale_y_continuous(breaks = pretty_breaks(n = 4)) +
    labs(
      title = NULL,   # column headings are added when the four panels are assembled
      x = paste0("PC1 (", round(var_exp[1], 1), "% variance)"),
      y = paste0("PC2 (", round(var_exp[2], 1), "% variance)"),
      color = NULL,
      shape = NULL
    ) +
    theme_minimal() +
    theme(
      legend.position = "top",
      legend.title = element_blank()
    )
}


# Matrices contain only analytes; metadata is exported alongside each stage.
# Pre and post metabolomics have different QC sample sets, so use their own metadata.
normalize_plate <- function(x) {
  x <- as.character(x)
  x[x == "S-trap batch 1"] <- "2"
  x[x == "S-trap batch 2"] <- "7"
  x <- ifelse(grepl("^[1-7]$", x), paste("Plate", x), x)
  gsub("^plate\\s*", "Plate ", x, ignore.case = TRUE)
}

pca_panel <- function(dataset, stage, title_text) {
  data_mat <- get(paste(stage, dataset, sep = "_"))
  metadata <- get(paste(stage, dataset, "metadata", sep = "_"))
  stopifnot(is.matrix(data_mat), is.numeric(data_mat),
            nrow(data_mat) == nrow(metadata),
            identical(rownames(data_mat), rownames(metadata)),
            all(is.finite(data_mat)))
  legend_group <- ifelse(metadata$sample_type == "SPQC", "SPQC",
                         normalize_plate(metadata$plate))
  stopifnot(!anyNA(legend_group), all(legend_group %in% legend_levels))
  make_pca_plot(data_mat, legend_group, paste(title_text, stage))
}

p.a26met <- pca_panel("met_a26", "post", "A26 Metabolomics")
p.a26pro <- pca_panel("pro_a26", "post", "A26 Proteomics")
p.a27met <- pca_panel("met_a27", "post", "A27 Metabolomics")
p.a27pro <- pca_panel("pro_a27", "post", "A27 Proteomics")

# -----------------------------
# Combine legend + panels
# -----------------------------


panels <- ggarrange(
  p.a26met + theme(legend.position = "none") + panel_margin,
  p.a26pro + theme(legend.position = "none") + panel_margin,
  p.a27met + theme(legend.position = "none") + panel_margin,
  p.a27pro + theme(legend.position = "none") + panel_margin,
  ncol = 2, nrow = 2,
  labels = c("a", "b", "c", "d"),
  font.label = list(size = 12, face = "bold"),
  label.x = 0, label.y = 1, hjust = 0, vjust = 1
)



legend_df <- data.frame(
  x = 1:length(legend_levels),
  y = 1,
  LegendGroup = factor(legend_levels, levels = legend_levels)
)

legend_plot <- ggplot(legend_df, aes(x, y)) +
  geom_point(
    aes(color = LegendGroup, shape = LegendGroup),
    size = 3,
    alpha = 0,          # hide points in the actual plot area
    show.legend = TRUE
  ) +
  scale_color_manual(
    values = legend_colors,
    breaks = legend_levels,
    limits = legend_levels,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = legend_shapes,
    breaks = legend_levels,
    limits = legend_levels,
    drop = FALSE
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      order = 1,
      override.aes = list(
        alpha = 1,
        size = pca_legend_pt,
        shape = unname(legend_shapes[legend_levels])
      )
    ),
    shape = "none"
  ) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = pca_legend_txt),
    legend.key.width = grid::unit(0.7, "cm"),
    legend.key.height = grid::unit(0.5, "cm"),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )



final_plot_post <- assemble_pca_figure(legend_plot, panels)





# Fig 2 (post-ComBat) — single TIFF per SR (was PNG + two broken vector calls
# referencing undefined `final_plot`; collapsed to one uniform TIFF export).
ggsave(
  file.path(output_dir, paste0("post_PCA_4panel_shared_legend.", fig_fmt)),
  plot = final_plot_post,
  width = fig_w_mm,
  height = fig_h_mm,
  units = "mm",
  dpi = fig_dpi,
  bg = "white",
  device = fig_fmt,
  compression = "lzw"
)





################################ Figures for Pre-Combat ################################################


p.a26met <- pca_panel("met_a26", "pre", "A26 Metabolomics")
p.a26pro <- pca_panel("pro_a26", "pre", "A26 Proteomics")
p.a27met <- pca_panel("met_a27", "pre", "A27 Metabolomics")
p.a27pro <- pca_panel("pro_a27", "pre", "A27 Proteomics")

# -----------------------------
# Combine legend + panels
# -----------------------------


panels <- ggarrange(
  p.a26met + theme(legend.position = "none") + panel_margin,
  p.a26pro + theme(legend.position = "none") + panel_margin,
  p.a27met + theme(legend.position = "none") + panel_margin,
  p.a27pro + theme(legend.position = "none") + panel_margin,
  ncol = 2, nrow = 2,
  labels = c("a", "b", "c", "d"),
  font.label = list(size = 12, face = "bold"),
  label.x = 0, label.y = 1, hjust = 0, vjust = 1
)



legend_df <- data.frame(
  x = 1:length(legend_levels),
  y = 1,
  LegendGroup = factor(legend_levels, levels = legend_levels)
)

legend_plot <- ggplot(legend_df, aes(x, y)) +
  geom_point(
    aes(color = LegendGroup, shape = LegendGroup),
    size = 3,
    alpha = 0,          # hide points in the actual plot area
    show.legend = TRUE
  ) +
  scale_color_manual(
    values = legend_colors,
    breaks = legend_levels,
    limits = legend_levels,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = legend_shapes,
    breaks = legend_levels,
    limits = legend_levels,
    drop = FALSE
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      order = 1,
      override.aes = list(
        alpha = 1,
        size = pca_legend_pt,
        shape = unname(legend_shapes[legend_levels])
      )
    ),
    shape = "none"
  ) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = pca_legend_txt),
    legend.key.width = grid::unit(0.7, "cm"),
    legend.key.height = grid::unit(0.5, "cm"),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )



final_plot_pre <- assemble_pca_figure(legend_plot, panels)





# S2 (pre-ComBat) — single TIFF, same parameters as Fig 2 for cross-figure consistency.
ggsave(
  file.path(output_dir, paste0("pre_PCA_4panel_shared_legend.", fig_fmt)),
  plot = final_plot_pre,
  width = fig_w_mm,
  height = fig_h_mm,
  units = "mm",
  dpi = fig_dpi,
  bg = "white",
  device = fig_fmt,
  compression = "lzw"
)




