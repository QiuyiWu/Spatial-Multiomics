# Run from the repository root: Rscript --vanilla Code/regenerate_Figure4.R [output_directory]
# Uses individual metabolites and the same loading-weighted protein clusters as Figure S5.
suppressPackageStartupMessages({library(here); library(cowplot); library(ggplot2)})
here::i_am('README.md')
args <- commandArgs(TRUE)
out <- if(length(args)) args[1] else here::here('Figure','Figure4')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
out <- normalizePath(out)
options(proteinfix.neighbor_figure_dir=out)
pdf(file.path(out,'diagnostics.pdf'))
run_animal <- function(animal) {
  e <- new.env(parent=globalenv())
  lines <- readLines(here::here('Code',animal,paste0('SampleSize_',animal,'_spearman.Rmd')))
  chunks <- list()
  for(start in grep('^```\\{r',lines)) {
    label <- sub('^```\\{r *','',lines[start]); label <- sub('[,}].*$','',label)
    end <- start + which(lines[(start+1):length(lines)]=='```')[1]
    chunks[[label]] <- lines[(start+1):(end-1)]
  }
  for(label in c('setup','lib','preparation',paste(tolower(animal),'met'),
                'scaling met','map met',paste(tolower(animal),'proteomic'),
                paste(tolower(animal),'pro'),'scaling','pro cluster')) {
    code <- chunks[[label]]
    # Only construct standardized inputs and loading-weighted clusters here.
    # Compute each full correlation matrix once below for the final maps.
    if(label %in% c('scaling met','scaling','pro cluster')) {
      cut <- grep('^cors\\.1 <-',code)[1]
      stopifnot(!is.na(cut));code <- code[seq_len(cut-1)]
    }
    message(animal,': ',label)
    eval(parse(text=code),e)
    if(label=='setup') {e$fig4_export <- FALSE; e$s5_export <- FALSE}
    if(label=='preparation') {
      A <- e$nbr1.full
      D <- igraph::distances(igraph::graph_from_adjacency_matrix(A,mode='undirected'))
      for(k in seq_along(e$neighbor_shells)) stopifnot(all(e$neighbor_shells[[k]]==(D==k)))
      message('Graph-distance validation passed: ',nrow(A),' samples')
    }
  }
  orders <- if(animal=='A26') 1:3 else c(1,3,5)
  correlations <- pair_correlations <- plots <- list()
  for(type in c('met','pro')) {
    X <- if(type=='met') e$X.sub.met else e$X.cluster.pro
    stopifnot(all(is.finite(X)),!anyDuplicated(rownames(X)))
    C <- cor(t(X),method='spearman')
    pair_correlations[[type]] <- list()
    plots[[type]] <- list()
    correlations[[type]] <- lapply(seq_along(orders),function(index) {
      k <- orders[index]
      A <- get(paste0('nbr',k,'.full'),e)[rownames(X),rownames(X),drop=FALSE]
      stopifnot(identical(rownames(X),rownames(A)),all(rowSums(A)>0))
      observed <- rowSums(C*A)/rowSums(A)
      # Independently evaluate complete neighbor means for representative samples.
      for(i in unique(round(seq(1,nrow(X),length.out=10)))) {
        neighbors <- which(A[i,]==1)
        direct <- mean(vapply(neighbors,function(j)
          cor(X[i,],X[j,],method='spearman'),numeric(1)))
        stopifnot(abs(observed[i]-direct)<1e-10)
      }
      pair_correlations[[type]][[index]] <<- C[which(A==1 & upper.tri(A),arr.ind=TRUE)]
      map <- e$paint_to_map(e$map.df,setNames(observed,rownames(X)))
      plots[[type]][[index]] <<- e$plot_map(map,NULL)
      observed
    })
  }
  result <- list(animal=animal,orders=orders,met=plots$met,pro=plots$pro,
                 correlations=correlations,pair_correlations=pair_correlations,
                 n=nrow(e$X.sub.met),met_features=ncol(e$X.sub.met),
                 protein_representation='loading_weighted_clusters',
                 protein_clusters=ncol(e$X.cluster.pro),
                 standardized_proteins=ncol(e$X.z.pro),constant_proteins=sum(!e$pro_keep))
  saveRDS(result[setdiff(names(result),c('met','pro'))],
          file.path(out,paste0(animal,'_Figure4_results.rds')))
  message(animal,': ',result$n,' spatial samples; ',result$met_features,
          ' metabolites; ',result$protein_clusters,' loading-weighted protein clusters')
  result
}
a26 <- run_animal('A26'); a27 <- run_animal('A27')
dev.off()
# Reference styling uses a shared blue-white-red correlation scale.
style_map <- function(p, order = NULL) {
  p + scale_fill_gradient2(low = "#4000FF", mid = "white", high = "#FF0000",
    midpoint = 0, limits = c(-1, 1), breaks = seq(-1, 1, .5),
    labels = sprintf("%.1f", seq(-1, 1, .5)), na.value = "grey90") +
    labs(title = order, fill = "avg cor") +
    theme(legend.position = "none", plot.title = element_text(hjust = 0, size = 12),
      axis.text = element_text(size = 3, colour = "grey55"),
      plot.margin = margin(2, 2, 2, 2))
}
make_group <- function(result,type,letter) {
  ordinals <- c("1st", "2nd", "3rd", "4th", "5th")
  plots <- lapply(seq_along(result$orders),function(i)
    style_map(result[[type]][[i]], ordinals[result$orders[i]]))
  triplet <- plot_grid(plotlist=plots,nrow=1)
  ggdraw() + draw_plot(triplet,x=.055,y=0,width=.945,height=.97) +
    draw_label(letter,x=.01,y=1,hjust=0,vjust=1,fontface='bold',size=18)
}
legend <- get_legend(style_map(a26$met[[1]]) +
 theme(legend.position='right',legend.title=element_text(size=9),
       legend.text=element_text(size=8), legend.key.width=grid::unit(3,'mm'),
       legend.key.height=grid::unit(5,'mm')))
body <- plot_grid(make_group(a26,'met','a'),make_group(a26,'pro','b'),
                  make_group(a27,'met','c'),make_group(a27,'pro','d'),ncol=2)
fig <- ggdraw() + draw_plot(body,x=.06,y=.075,width=.865,height=.865) +
 draw_label('Metabolomics',x=.28,y=.97,size=16) +
 draw_label('Protein clusters',x=.71,y=.97,size=16) +
 draw_label('A26',x=.035,y=.73,size=15) +
 draw_label('A27',x=.035,y=.295,size=15) +
 draw_label('Neighbor order',x=.5,y=.035,size=16) +
 draw_plot(legend,x=.925,y=.34,width=.075,height=.35)
ggsave(file.path(out,'Figure4_corrected.tiff'),fig,width=15,height=8,units='in',dpi=600,compression='lzw',bg='white')
ggsave(file.path(out,'Figure4_corrected.png'),fig,width=15,height=8,units='in',dpi=180,bg='white')
