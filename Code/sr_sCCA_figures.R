#------------------- A26 -----------------------------


srcca_component <- function(i, scca_res, Xtc, Ytc, map.df, X, Y, a26.pro.cluster, a26.met.cluster) {
  ai <- scca_res$u[, i]
  bi <- scca_res$v[, i]

  u_scores <- Xtc %*% ai
  v_scores <- Ytc %*% bi
  names(u_scores) <- rownames(X)
  names(v_scores) <- rownames(Y)

  u_unsmoothed <- W_half_inv %*% u_scores
  v_unsmoothed <- W_half_inv %*% v_scores
  names(u_unsmoothed) <- rownames(X)
  names(v_unsmoothed) <- rownames(Y)


  # Spatial maps
  u_matrix <- matrix(u_scores[match(map.df, names(u_scores))], nrow = nrow(map.df))
  u_df <- reshape2::melt(apply(u_matrix, 2, rev), varnames = c("Row", "Column"),
                         value.name = paste0("spCCA_U", i))
  v_matrix <- matrix(v_scores[match(map.df, names(v_scores))], nrow = nrow(map.df))
  v_df <- reshape2::melt(apply(v_matrix, 2, rev), varnames = c("Row", "Column"),
                         value.name = paste0("spCCA_V", i))

  hjust_u <- ifelse(i < 10, -7, -4.5)
  hjust_v <- ifelse(i < 10, -7, -4.5)
  vjust_u <- ifelse(i < 10,  1.5,  1.5)
  vjust_v <- ifelse(i < 10,  1.5,  1.5)

  p_u <- ggplot(u_df, aes_string("Column", "Row", fill = paste0("spCCA_U", i))) +
    geom_raster() +
    scale_x_reverse() +
    coord_equal() +
    scale_fill_viridis_c(option = "plasma", na.value = "grey90") +
    annotate(
      "text",
      x = max(u_df$Column, na.rm = TRUE),
      y = max(u_df$Row, na.rm = TRUE),
      label = paste0("U", i),
      hjust = hjust_u, vjust = vjust_u,
      size = 8
    ) +
    labs(x = NULL, y = NULL, fill = "Score") +
    theme_minimal() +
    theme(panel.grid = element_blank())

  p_v <- ggplot(v_df, aes_string("Column", "Row", fill = paste0("spCCA_V", i))) +
    geom_raster() +
    scale_x_reverse() +
    coord_equal() +
    scale_fill_viridis_c(option = "plasma", na.value = "grey90") +
    annotate(
      "text",
      x = max(u_df$Column, na.rm = TRUE),
      y = max(u_df$Row, na.rm = TRUE),
      label = paste0("V", i),
      hjust = hjust_v, vjust = vjust_v,
      size = 8
    ) +
    labs(x = NULL, y = NULL, fill = "Score") +
    theme_minimal() +
    theme(panel.grid = element_blank())

  print(p_v + p_u)
}




srcca_component(2, scca_res, Xtc, Ytc, map.df, X, Y, a26.pro.cluster, a26.met.cluster)



dir.create("A26_components", showWarnings = FALSE)
for (i in 1:10) {
  final_plot_pre <- srcca_component(
    i, scca_res, Xtc, Ytc, map.df, X, Y,
    a26.pro.cluster, a26.met.cluster
  )

  ggsave(
    filename = file.path("A26_components", paste0("A26_comp", i, ".png")),
    plot = final_plot_pre,
    width = 25,
    height = 20,
    units = "cm",
    dpi = 300
  )
}















#------------------- A27 -----------------------------


srcca_component <- function(i, scca_res, Xtc, Ytc, map.df, X, Y, a27.pro.cluster, a27.met.cluster) {
  ai <- scca_res$u[, i]
  bi <- scca_res$v[, i]

  u_scores <- Xtc %*% ai
  v_scores <- Ytc %*% bi
  names(u_scores) <- rownames(X)
  names(v_scores) <- rownames(Y)

  u_unsmoothed <- W_half_inv %*% u_scores
  v_unsmoothed <- W_half_inv %*% v_scores
  names(u_unsmoothed) <- rownames(X)
  names(v_unsmoothed) <- rownames(Y)


  # Spatial maps
  u_matrix <- matrix(u_scores[match(map.df, names(u_scores))], nrow = nrow(map.df))
  u_df <- reshape2::melt(apply(u_matrix, 2, rev), varnames = c("Row", "Column"),
                         value.name = paste0("spCCA_U", i))
  v_matrix <- matrix(v_scores[match(map.df, names(v_scores))], nrow = nrow(map.df))
  v_df <- reshape2::melt(apply(v_matrix, 2, rev), varnames = c("Row", "Column"),
                         value.name = paste0("spCCA_V", i))

  p_u <- ggplot(u_df, aes_string("Column", "Row", fill = paste0("spCCA_U", i))) +
    geom_raster() +
    coord_equal() +
    scale_fill_viridis_c(option = "plasma", na.value = "grey90") +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = paste0("U", i),
      hjust = 1.5, vjust = 2.5,
      size = 8
    ) +
    labs(x = NULL, y = NULL, fill = "Score") +
    theme_minimal() +
    theme(panel.grid = element_blank())

  p_v <- ggplot(v_df, aes_string("Column", "Row", fill = paste0("spCCA_V", i))) +
    geom_raster() +
    coord_equal() +
    scale_fill_viridis_c(option = "plasma", na.value = "grey90") +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = paste0("V", i),
      hjust = 1.5, vjust = 2.5,
      size = 8
    ) +
    labs(x = NULL, y = NULL, fill = "Score") +
    theme_minimal() +
    theme(panel.grid = element_blank())

  print(p_v + p_u)
  }




srcca_component(1, scca_res, Xtc, Ytc, map.df, X, Y, a27.pro.cluster, a27.met.cluster)
dir.create("A27_components", showWarnings = FALSE)
for (i in 1:10) {
  final_plot_pre <- srcca_component(
    i, scca_res, Xtc, Ytc, map.df, X, Y,
    a27.pro.cluster, a27.met.cluster
  )

  ggsave(
    filename = file.path("A27_components", paste0("A27_comp", i, ".png")),
    plot = final_plot_pre,
    width = 25,
    height = 20,
    units = "cm",
    dpi = 300
  )
}
