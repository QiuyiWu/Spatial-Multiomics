# Exact row removal: sample labels can repeat between animals.
figure5_outliers <- function(data, scores) {
  stopifnot(nrow(data) == nrow(scores),
            identical(as.character(data$Sample), as.character(scores$Sample)),
            identical(as.character(data$Subject), as.character(scores$Subject)),
            is.logical(scores$outlier), !anyNA(scores$outlier))
  key <- paste(data$Subject, data$Sample, sep = '::')
  if(anyDuplicated(key)) stop('Duplicate Subject/Sample key in Figure 5 biological samples')
  legacy_removed <- data$Sample %in% scores$Sample[scores$outlier]
  data.frame(Subject = data$Subject, Sample = data$Sample,
             flagged = scores$outlier, legacy_removed = legacy_removed,
             extra_legacy_removal = legacy_removed & !scores$outlier)
}

# Normalize both sides without losing original labels; reject ambiguous keys.
figure5_match_loadings <- function(pc1, weights) {
  pc_keys <- make.names(pc1$Feature)
  weight_keys <- make.names(rownames(weights))
  if(anyDuplicated(pc_keys) || anyDuplicated(weight_keys))
    stop('Ambiguous metabolite names after make.names normalization')
  common <- intersect(pc_keys, weight_keys)
  if(length(common) < 3L) stop('Too few shared metabolites for loading comparisons')
  p <- pc1[match(common, pc_keys), , drop = FALSE]
  p$OriginalFeature <- p$Feature
  p$Feature <- common
  rownames(p) <- common
  w <- weights[match(common, weight_keys), , drop = FALSE]
  rownames(w) <- common
  stopifnot(identical(p$Feature, rownames(w)))
  list(pc1 = p, weights = w,
       mapping = data.frame(key = common, PCA_feature = p$OriginalFeature,
                            cluster_feature = rownames(weights)[match(common, weight_keys)]),
       unmatched = pc1$Feature[!pc_keys %in% common],
       literal_matches = length(intersect(pc1$Feature, rownames(weights))))
}
