# Use current C5 loadings rather than a separately exported contribution workbook.
protein_cluster_contributions <- function(animal) {
  e <- new.env(parent=emptyenv())
  load(here::here('Data',paste0(animal,'_PChclust.RData')),envir=e)
  L <- get(paste0(tolower(animal),'.pro.cluster'),e)$loading
  idx <- which(L!=0,arr.ind=TRUE)
  data.frame(Feature=rownames(L)[idx[,1]],Cluster=colnames(L)[idx[,2]],
             Loading=L[idx],stringsAsFactors=FALSE)
}
