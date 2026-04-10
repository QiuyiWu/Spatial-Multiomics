library(ggplot2)
library(ggpubr)
library(scales)

load("data/met_a26_PCAdata.RData")
load("data/met_a27_PCAdata.RData")
load("data/pro_a26_PCAdata.RData")
load("data/pro_a27_PCAdata.RData")

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

################################ Figures for Post-Combat ################################################

# -----------------------------
# Common plotting function
# -----------------------------
make_pca_plot <- function(data_mat, legend_group, title_text) {
  pca_res <- prcomp(data_mat, center = TRUE, scale. = TRUE)
  pca_df <- as.data.frame(pca_res$x[, 1:2])
  colnames(pca_df) <- c("PC1", "PC2")
  
  var_exp <- 100 * (pca_res$sdev^2 / sum(pca_res$sdev^2))
  
  pca_df$LegendGroup <- factor(legend_group, levels = legend_levels)
  
  ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = LegendGroup, shape = LegendGroup), size = 2) +
    scale_color_manual(values = legend_colors, drop = FALSE) +
    scale_shape_manual(values = legend_shapes, drop = FALSE) +
    scale_x_continuous(breaks = pretty_breaks(n = 4)) +
    scale_y_continuous(breaks = pretty_breaks(n = 4)) +
    labs(
      title = title_text,
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


pre_met_a26 <- pre_met_a26[-c(404:406), ]

plate_a26met <- pre_met_a26$color_indicator
plate_a26met <- gsub("^plate", "Plate", plate_a26met)
plate_a26met <- gsub("^Plate\\s*", "Plate ", plate_a26met)

legend_a26met <- ifelse(pre_met_a26$sample_type == "SPQC", "SPQC", plate_a26met)

p.a26met <- make_pca_plot(
  data_mat = post_met_a26,
  legend_group = legend_a26met,
  title_text = "A26 Metabolomics"
)

plate_a26pro <- pre_pro_a26$pro.plate
plate_a26pro <- gsub("^Plate", "Plate ", plate_a26pro)   # Plate1 -> Plate 1

legend_a26pro <- ifelse(pre_pro_a26$pro.sampletype == "SPQC", "SPQC", plate_a26pro)

p.a26pro <- make_pca_plot(
  data_mat = post_pro_a26,
  legend_group = legend_a26pro,
  title_text = "A26 Proteomics"
)

plate_a27met <- pre_met_a27$color_indicator
plate_a27met <- gsub("^plate", "Plate", plate_a27met)
plate_a27met <- gsub("^Plate\\s*", "Plate ", plate_a27met)

type_a27met <- ifelse(grepl("SPQC", pre_met_a27$sampleID), "SPQC", "Non-SPQC")
legend_a27met <- ifelse(type_a27met == "SPQC", "SPQC", plate_a27met)

p.a27met <- make_pca_plot(
  data_mat = post_met_a27,
  legend_group = legend_a27met,
  title_text = "A27 Metabolomics"
)

plate_a27pro <- as.character(pre_pro_a27$pro.plate)

plate_a27pro[plate_a27pro == "1"] <- "Plate 1"
plate_a27pro[plate_a27pro == "3"] <- "Plate 3"
plate_a27pro[plate_a27pro == "4"] <- "Plate 4"
plate_a27pro[plate_a27pro == "5"] <- "Plate 5"
plate_a27pro[plate_a27pro == "6"] <- "Plate 6"
plate_a27pro[plate_a27pro == "S-trap batch 1"] <- "Plate 2"
plate_a27pro[plate_a27pro == "S-trap batch 2"] <- "Plate 7"

legend_a27pro <- ifelse(pre_pro_a27$pro.sampletype == "SPQC", "SPQC", plate_a27pro)

p.a27pro <- make_pca_plot(
  data_mat = post_pro_a27,
  legend_group = legend_a27pro,
  title_text = "A27 Proteomics"
)


# -----------------------------
# Combine legend + panels
# -----------------------------


panels <- ggarrange(
  p.a26met + theme(legend.position = "none"),
  p.a26pro + theme(legend.position = "none"),
  p.a27met + theme(legend.position = "none"),
  p.a27pro + theme(legend.position = "none"),
  ncol = 2, nrow = 2,
  labels = c("A", "B", "C", "D")
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
        size = 3,
        shape = unname(legend_shapes[legend_levels])
      )
    ),
    shape = "none"
  ) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.key.width = grid::unit(0.7, "cm"),
    legend.key.height = grid::unit(0.5, "cm"),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )



final_plot_post <- ggarrange(
  legend_plot,
  panels,
  ncol = 1,
  heights = c(0.08, 1)
)





ggsave(
  "post_PCA_4panel_shared_legend.png",
  plot = final_plot_post,
  width = 25,
  height = 20,
  units = "cm",
  dpi = 300
)



ggsave(
  "PCA_4panel_shared_legend.pdf",
  plot = final_plot,
  width = 25,
  height = 20,
  units = "cm"
)

ggsave(
  "PCA_4panel_shared_legend.eps",
  plot = final_plot,
  width = 25,
  height = 20,
  units = "cm",
  device = cairo_ps
)





################################ Figures for Pre-Combat ################################################


p.a26met <- make_pca_plot(
  data_mat = pre_met_a26[, -c(1:3)],
  legend_group = legend_a26met,
  title_text = "A26 Metabolomics"
)


p.a26pro <- make_pca_plot(
  data_mat = pre_pro_a26[, -c(1:3)],
  legend_group = legend_a26pro,
  title_text = "A26 Proteomics"
)


p.a27met <- make_pca_plot(
  data_mat = pre_met_a27[, -c(631:633)],
  legend_group = legend_a27met,
  title_text = "A27 Metabolomics"
)

p.a27pro <- make_pca_plot(
  data_mat = pre_pro_a27[, -c(1:4)],
  legend_group = legend_a27pro,
  title_text = "A27 Proteomics"
)


# -----------------------------
# Combine legend + panels
# -----------------------------


panels <- ggarrange(
  p.a26met + theme(legend.position = "none"),
  p.a26pro + theme(legend.position = "none"),
  p.a27met + theme(legend.position = "none"),
  p.a27pro + theme(legend.position = "none"),
  ncol = 2, nrow = 2,
  labels = c("A", "B", "C", "D")
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
        size = 3,
        shape = unname(legend_shapes[legend_levels])
      )
    ),
    shape = "none"
  ) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.key.width = grid::unit(0.7, "cm"),
    legend.key.height = grid::unit(0.5, "cm"),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )



final_plot_pre <- ggarrange(
  legend_plot,
  panels,
  ncol = 1,
  heights = c(0.08, 1)
)





ggsave(
  "pre_PCA_4panel_shared_legend.png",
  plot = final_plot_pre,
  width = 25,
  height = 20,
  units = "cm",
  dpi = 300
)




