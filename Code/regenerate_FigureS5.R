# Run from the repository root: Rscript --vanilla Code/regenerate_FigureS5.R [output_directory]
suppressPackageStartupMessages(library(here))
here::i_am('README.md')
args <- commandArgs(TRUE)
out <- if(length(args)) args[1] else here::here('Figure','FigureS5')
dir.create(out,recursive=TRUE,showWarnings=FALSE);out<-normalizePath(out)
options(proteinfix.neighbor_figure_dir=out)
pdf(file.path(out,'S5_diagnostics.pdf'))
run <- function(animal) {
 e <- new.env(parent=globalenv())
 lines <- readLines(here::here('Code',animal,paste0('SampleSize_',animal,'_spearman.Rmd')))
 chunks<-list()
 for(start in grep('^```\\{r',lines)) {
  label<-sub('^```\\{r *','',lines[start]);label<-sub('[,}].*$','',label)
  end<-start+which(lines[(start+1):length(lines)]=='```')[1]
  chunks[[label]]<-lines[(start+1):(end-1)]
 }
 for(label in c('setup','lib','preparation',paste(tolower(animal),'met'),'scaling met',
                paste(tolower(animal),'proteomic'),paste(tolower(animal),'pro'),'scaling','pro cluster')) {
  code<-chunks[[label]]
  # The S5 runner computes pair correlations once rather than repeating per-node calls.
  if(label %in% c('scaling met','scaling','pro cluster')) {
   cut<-grep('^cors\\.1 <-',code)[1]
   code<-code[seq_len(cut-1)]
  }
  message(animal,': ',label);eval(parse(text=code),e)
  if(label=='setup') {e$fig4_export<-FALSE;e$s5_export<-FALSE}
 }
 orders<-if(animal=='A26') 1:3 else c(1,3,5)
 result<-list(animal=animal,orders=orders,n=nrow(e$X.sub.met),clusters=ncol(e$X.cluster.pro))
 for(type in c('met','procluster')) {
  X<-if(type=='met') e$X.sub.met else e$X.cluster.pro
  stopifnot(all(is.finite(X)),ncol(X)>1)
  C<-cor(t(X),method='spearman')
  result[[type]]<-lapply(orders,function(k) {
   A<-get(paste0('nbr',k,'.full'),e)[rownames(X),rownames(X),drop=FALSE]
   idx<-which(A==1 & upper.tri(A),arr.ind=TRUE)
   values<-C[idx];stopifnot(all(is.finite(values)))
   # Independent direct-pair check on a representative subset.
   for(j in unique(round(seq(1,nrow(idx),length.out=10))))
    stopifnot(abs(values[j]-cor(X[idx[j,1],],X[idx[j,2],],method='spearman'))<1e-12)
   values
  })
 }
 saveRDS(result,file.path(out,paste0(animal,'_FigureS5_results.rds')))
 message(animal,': ',result$n,' samples; ',result$clusters,' protein clusters')
 result
}
a26<-run('A26');a27<-run('A27');dev.off()
# Base-R boxes: median bars, 1.5-IQR whiskers, no outlier symbols; dots/lines are means.
render <- function() {
 par(mfrow=c(2,2),mar=c(3,3,1.5,1),oma=c(4,5,3,1),mgp=c(1.8,.6,0),tcl=-.3,
     cex.axis=.8,bty='o',col.axis='grey30',col='grey60')
 ordinal<-c('1st','2nd','3rd','4th','5th')
 for(r in list(a26,a27)) for(type in c('met','procluster')) {
  x<-r[[type]];ylim<-range(unlist(x),finite=TRUE)
  boxplot(x,names=ordinal[r$orders],outline=FALSE,col='grey88',border='grey65',
          boxwex=.45,ylim=ylim,ylab='',xlab='',whisklty=2,medcol='grey40',medlwd=1)
  lines(seq_along(x),vapply(x,mean,numeric(1)),type='b',pch=16,cex=.5,col='grey35',lwd=.8)
 }
 mtext('Metabolomic order',side=3,outer=TRUE,at=.25,line=1,font=2,cex=1.1,col='black')
 mtext('Proteomic cluster similarity',side=3,outer=TRUE,at=.75,line=1,font=2,cex=1.1,col='black')
 mtext('A26',side=2,outer=TRUE,at=.77,line=1,font=2,las=1,cex=1,col='black')
 mtext('A27',side=2,outer=TRUE,at=.27,line=1,font=2,las=1,cex=1,col='black')
 mtext('Average correlation',side=2,outer=TRUE,line=3.5,font=2,cex=1.2,col='black')
 mtext('Neighbor order',side=1,outer=TRUE,line=2,font=2,cex=1.2,col='black')
}
tiff(file.path(out,'FigureS5_corrected.tiff'),width=10,height=6.4,units='in',res=600,type='cairo',compression='lzw');render();dev.off()
png(file.path(out,'FigureS5_corrected.png'),width=10,height=6.4,units='in',res=180);render();dev.off()
