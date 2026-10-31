# Run from the repository root: Rscript --vanilla Code/regenerate_Figure4.R [output_directory]
# Executes the raw metabolite/protein Figure 4 chunks, excluding the optional cluster analysis.
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
  starts <- grep('^```\\{r',lines)
  for(start in starts) {
    label <- sub('^```\\{r *','',lines[start]); label <- sub('[,}].*$','',label)
    if(label=='pro cluster') break
    end <- start + which(lines[(start+1):length(lines)]=='```')[1]
    message(animal,': ',label)
    eval(parse(text=lines[(start+1):(end-1)]),e)
    if(label=='setup') {e$fig4_export <- FALSE; e$s5_export <- FALSE}
    if(label=='preparation') {
      A <- e$nbr1.full
      D <- igraph::distances(igraph::graph_from_adjacency_matrix(A,mode='undirected'))
      for(k in seq_along(e$neighbor_shells)) stopifnot(all(e$neighbor_shells[[k]]==(D==k)))
      message('Graph-distance validation passed: ',nrow(A),' samples')
    }
  }
  orders <- if(animal=='A26') 1:3 else c(1,3,5)
  correlations <- list()
  for(type in c('met','pro')) {
    X <- if(type=='met') e$X.sub.met else e$X.z.pro
    stopifnot(all(is.finite(X)))
    C <- cor(t(X),method='spearman')
    correlations[[type]] <- lapply(orders,function(k) {
      A <- get(paste0('nbr',k,'.full'),e)[rownames(X),rownames(X),drop=FALSE]
      observed <- get(paste0('cor',k,'.',type),e)
      expected <- rowSums(C*A)/rowSums(A)
      stopifnot(identical(rownames(X),rownames(A)),max(abs(observed-expected),na.rm=TRUE)<1e-10)
      observed
    })
  }
  result <- list(animal=animal,orders=orders,met=get(paste0('p_',tolower(animal),'_met'),e),
                 pro=get(paste0('p_',tolower(animal),'_pro'),e),correlations=correlations,
                 n=nrow(e$X.sub.met),met_features=ncol(e$X.sub.met),
                 protein_features=ncol(e$X.z.pro),constant_proteins=sum(!e$pro_keep))
  saveRDS(result[setdiff(names(result),c('met','pro'))],
          file.path(out,paste0(animal,'_Figure4_results.rds')))
  message(animal,': ',result$n,' spatial samples; ',result$met_features,' metabolites; ',result$protein_features,' variable proteins; ',result$constant_proteins,' constant proteins omitted')
  result
}
a26 <- run_animal('A26'); a27 <- run_animal('A27')
dev.off()
make_group <- function(result,type,letter) {
  plots <- lapply(seq_along(result$orders),function(i) result[[type]][[i]] +
    labs(title=as.character(result$orders[i])) +
    theme(legend.position='none',plot.title=element_text(hjust=.5,size=12),plot.margin=margin(2,2,2,2)))
  triplet <- plot_grid(plotlist=plots,nrow=1)
  ggdraw() + draw_plot(triplet,x=.025,y=0,width=.975,height=.88) +
    draw_label(letter,x=.01,y=.99,hjust=0,vjust=1,fontface='bold',size=18) +
    draw_label('Neighbor order',x=.51,y=.96,size=11)
}
legend <- get_legend(a26$met[[1]] + labs(fill='Average\nSpearman\ncorrelation') +
 theme(legend.position='right',legend.title=element_text(size=10),legend.text=element_text(size=9)))
body <- plot_grid(make_group(a26,'met','a'),make_group(a26,'pro','b'),
                  make_group(a27,'met','c'),make_group(a27,'pro','d'),ncol=2)
fig <- ggdraw() + draw_plot(body,x=.055,y=.04,width=.865,height=.89) +
 draw_label('Metabolomics',x=.27,y=.97,size=16) +
 draw_label('Proteomics',x=.70,y=.97,size=16) +
 draw_label('A26',x=.025,y=.705,size=15,fontface='bold',angle=90) +
 draw_label('A27',x=.025,y=.26,size=15,fontface='bold',angle=90) +
 draw_plot(legend,x=.92,y=.32,width=.08,height=.35)
ggsave(file.path(out,'Figure4_corrected.tiff'),fig,width=15,height=9,units='in',dpi=600,compression='lzw',bg='white')
ggsave(file.path(out,'Figure4_corrected.png'),fig,width=15,height=9,units='in',dpi=180,bg='white')
