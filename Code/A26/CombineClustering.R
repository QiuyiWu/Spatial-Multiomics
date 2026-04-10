
library(Hmisc)
library(sRDA)
library(PMA)
library(corrplot)
library(dplyr)
library(readxl)
library(ggplot2)
library(ggtext)
library(pheatmap)
library(gridExtra)
library(factoextra)
library(kableExtra)
library(reshape2)
library(pROC)
library(randomForest)
library(glmnet)
library(DT)
library(reshape2)
library(sva)
library(writexl)
library(ggrepel)
library(purrr)
library(tidyverse)
library(plyr)
library(DT)
library(htmltools)

##################################################
# C4. Selecting Mapping Samples ######
##################################################
# finally 267 samples
a26.met.data <- as.data.frame(read_excel("../data/PostCombat_A26_Metabolite_beforeHeldout.xlsx"))
a26.pro.data <- as.data.frame(read_excel("../data/PostCombat_A26_Protein.xlsx", sheet = "combat_corrected_unscaled"))
# convert Duke ID to Analysis ID
dukeid.df <- data.frame(read_excel("../data/A26_Bonnie_Flatmap_Grid_20260315.xlsx", sheet = "Proteome_DukeIDs"))[, -1]
colnames(a26.pro.data)[1] <- "DukeID"
a26.pro.data$DukeID <- as.numeric(substr(a26.pro.data$DukeID, 3,8))
# Join and replace DukeID with AnalysisID
a26.pro.data <- a26.pro.data %>%
  left_join(dukeid.df, by = "DukeID") %>%
  mutate(DukeID = AnalysisID) %>%       
  select(-AnalysisID)                   

heldout_updated <- as.data.frame(read_excel("../data/Heldout_Updated_a26.xlsx"))

a26.met.data <- a26.met.data[a26.met.data$Sample %in% unique(heldout_updated$finalID ), ]
a26.pro.data <- a26.pro.data[a26.pro.data$DukeID %in% unique(heldout_updated$finalID ), ]

rownames(a26.met.data) <- a26.met.data$Sample
rownames(a26.pro.data) <- a26.pro.data$DukeID
# save(a26.met.data, a26.pro.data, file = "a26.RData")


##################################################
# C5. Clustering Analysis
##################################################

# recursive splitting by both variance explained by PC1 and correlation between two clusters
PChclust <- function(x, method.cor = "spearman", method.hclust="average",
                     min_pc1_var = 0.75, max_pc1_cor = 0.7, min_cluster_size = 1) {
  # Scale variables
  x <- scale(x)
  
  # Compute variable correlations & hierarchical clustering tree
  if (method.hclust=="hclustvar"){
    clustvar_res <- hclustvar(X.quanti = x)
    dend <- as.dendrogram(as.hclust(clustvar_res))
  }else{
    corx <- cor(x, method=method.cor)
    distx <- as.dist(1 - corx^2)
    treex <- hclust(distx, method=method.hclust)
    dend <- as.dendrogram(treex)
  }
  
  # Output group vector
  group <- rep("cls", ncol(x))
  names(group) <- colnames(x)
  
  # Recursive helper
  # Pre-calculate column indices to avoid string matching labels
  col_indices <- setNames(seq_len(ncol(x)), colnames(x))
  
  assign_groups <- function(subdend, prefix) {
    lbls <- labels(subdend)
    idx <- col_indices[lbls]
    
    s_values <- svd(x[, idx, drop = FALSE], nu = 0, nv = 0)$d
    var_explained_pc1 <- (s_values[1]^2) / sum(s_values^2)
    
    if (var_explained_pc1 >= min_pc1_var || is.leaf(subdend) || length(lbls) <= min_cluster_size) {
      group[lbls] <<- prefix
    } else {
      # Recursion continues, but the memory 'baggage' is now just a few integers
      assign_groups(subdend[[1]], paste0(prefix, ".1"))
      assign_groups(subdend[[2]], paste0(prefix, ".2"))
    }
  }
  
  assign_groups(dend, "cls")
  return(group)
}


PChclust.score <- function(x, group){
  grp <- table(group)
  ngrp <- length(grp)
  # rename clusters
  name_grp <- names(grp)
  for (ii in 1:ngrp){
    group[group==name_grp[ii]] <- paste0("Cluster",ii)
  }
  grp <- table(group)
  
  score <- matrix(0, nrow(x), ngrp)
  rownames(score) <- rownames(x)
  loading <- matrix(0, ncol(x), ngrp)
  rownames(loading) <- colnames(x)
  colnames(score) <- colnames(loading) <- names(grp)
  var_pct <- rep(0,ngrp)
  names(var_pct) <- names(grp)
  
  # get the loading and score of the top PC in each cluster
  for (ii in 1:length(grp)){
    var_ii <- group==names(grp)[ii]
    x_ii <- x[,var_ii]
    res_svd <- svd(scale(x_ii), nu=1, nv=1)
    score[,ii] <- res_svd$u[,1]
    loading[var_ii,ii] <- res_svd$v[,1]
    var_pct[ii] <- res_svd$d[1]^2/sum(res_svd$d^2)
  }
  
  # reorder the loading and scores by variable cluster
  loading_chunks <- vector("list", ncol(loading))
  for(ii in seq_len(ncol(loading))) {
    tmp <- loading[loading[, ii] != 0, , drop = FALSE]
    tmp <- tmp[order(tmp[, ii]), , drop = FALSE]
    loading_chunks[[ii]] <- tmp
  }
  loading_reorder <- do.call(rbind, loading_chunks)
  
  # reorder the score accordingly
  score <- score[,colnames(loading)]
  for(ii in 1:ncol(loading)){
    tmp <- loading[loading[,ii]!=0,,drop=FALSE]
    if(nrow(tmp)==1){
      colnames(loading_reorder)[ii] <- rownames(tmp)
      colnames(score)[ii] <- rownames(tmp)
    }
  }
  num_var <- as.data.frame(table(group))
  colnames(num_var) <- c("Cluster", "num_var")
  num_var$var_pct <- var_pct[num_var[,1]]
  return(list(group=group, score=score, loading=loading_reorder, clust_summary=num_var, var_pct=var_pct))
}


### LOAD DATA
met_data <- a26.met.data 
pro_data <- a26.pro.data 


# MET data
met_matrix <- data.frame(met_data)
rownames(met_matrix) <- met_data$Sample
met_matrix <- met_matrix[,-(1:3)]
met_matrix <- scale(met_matrix)
x <- met_matrix
group <- PChclust(x, method.hclust="average", min_pc1_var=0.7, min_cluster_size = 1)
length(table(group))
range(table(group))
a26.met.cluster <- PChclust.score(x, group)

# PRO data
pro_matrix <- data.frame(pro_data)
rownames(pro_matrix) <- pro_matrix$AnalysisID
pro_matrix <- pro_matrix[,-(1:3)]
# remove constant variable
pro_matrix <- pro_matrix[,apply(pro_matrix,2,sd)!=0]
pro_matrix <- scale(pro_matrix)
x <- pro_matrix
group <- PChclust(x, method.hclust="average", min_pc1_var=0.7, min_cluster_size = 1)
a26.pro.cluster <- PChclust.score(x, group)

#save(a26.met.cluster, a26.pro.cluster, file = "../data/A26_PChclust.RData")




