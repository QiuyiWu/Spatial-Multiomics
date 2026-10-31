# Run from the repository root, for example:
# Rscript --vanilla Code/regenerate_C6_outputs.R A26 /path/to/CorrectedFigure
# Recomputes the spatial C6 model (seed 123, 25 permutations, 10 components),
# then runs the corresponding section of the supplied sr_sCCA_figures.R.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2)
animal <- toupper(args[[1]])
stopifnot(animal %in% c('A26','A27'))
output_root <- normalizePath(args[[2]], mustWork = FALSE)
repo <- if(length(args)>=3) normalizePath(args[[3]]) else normalizePath(getwd())
setwd(repo)
here::i_am('README.md')
dir.create(output_root, recursive=TRUE, showWarnings=FALSE)
diagnostic_dir <- here::here('Figure','C6_diagnostics')
dir.create(diagnostic_dir,recursive=TRUE,showWarnings=FALSE)
pdf(file.path(diagnostic_dir,paste0(animal,'_C6_diagnostics.pdf')))
lines <- readLines(file.path('Code',animal,paste0(animal,'_analysis_sCCA_laplacians.Rmd')))
starts <- grep('^```\\{r',lines)
chunks <- list()
for(i in starts) {
  label <- sub('^```\\{r[[:space:]]*','',lines[i])
  label <- trimws(sub('[,}].*$','',label))
  end <- i + which(lines[(i+1):length(lines)]=='```')[[1]]
  chunks[[label]] <- lines[(i+1):(end-1)]
}
graph_chunk <- 'spCCA cluster 1st component'
for(label in c('setup','omics map',graph_chunk,'sparse spCCA 2')) {
  stopifnot(label %in% names(chunks))
  message(animal,': running ',label,' at ',format(Sys.time()))
  eval(parse(text=chunks[[label]]),envir=.GlobalEnv)
  if(label=='omics map') {
    stopifnot(!anyNA(xy_ordered),!anyDuplicated(rownames(xy_ordered)),
              !anyDuplicated(xy_ordered))
    # Maps can show one sample in several cells. Preserve the Rmds' match()
    # rule (first coordinate per sample) for the graph, and render all cells.
    message(animal,': ',sum(duplicated(xy_coords$sample)),
            ' repeated map cells; graph uses one coordinate per sample')
    message(animal,': ',nrow(xy_ordered),' spatially matched samples')
  }
}
stopifnot(ncol(scca_res$u)==10,ncol(scca_res$v)==10,
          nrow(scca_res$u)==ncol(Xtc),nrow(scca_res$v)==ncol(Ytc),
          identical(rownames(X),rownames(Y)),identical(rownames(X),rownames(xy_ordered)),
          all(is.finite(scca_res$u)),all(is.finite(scca_res$v)),
          all(is.finite(W_half)),all(is.finite(W_half_inv)),
          max(abs(W_half %*% W_half_inv-diag(nrow(W_half))))<1e-8,
          isSymmetric(W_raw),all(diag(W_raw)==0),all(rowSums(W_raw)>0))
u_scores <- Xtc %*% scca_res$u
v_scores <- Ytc %*% scca_res$v
stopifnot(all(is.finite(u_scores)),all(is.finite(v_scores)))
component_correlations <- vapply(1:10,function(i)cor(u_scores[,i],v_scores[,i]),numeric(1))
stopifnot(all(is.finite(component_correlations)))
summary <- data.frame(animal=animal,component=1:10,samples=nrow(X),
                      protein_clusters=ncol(X),metabolite_clusters=ncol(Y),
                      penaltyx=penaltyx,penaltyz=penaltyz,
                      map_score_correlation=component_correlations)
write.csv(summary,file.path(diagnostic_dir,paste0(animal,'_C6_summary.csv')),row.names=FALSE)
saveRDS(list(scca_res=scca_res,Xtc=Xtc,Ytc=Ytc,X=X,Y=Y,map.df=map.df,
             xy_ordered=xy_ordered,W_half=W_half,W_half_inv=W_half_inv,
             W_raw=W_raw,perm_out=perm_out,penaltyx=penaltyx,penaltyz=penaltyz),
        file.path(diagnostic_dir,paste0(animal,'_C6_state.rds')))
# Execute the supplied plotting code in the freshly computed model workspace.
figures <- readLines(here::here('Code','sr_sCCA_figures.R'))
separator <- grep('^#-+ A27 ',figures)
stopifnot(length(separator)==1)
figure_lines <- if(animal=='A26')figures[seq_len(separator-1)] else figures[separator:length(figures)]
figure_dir <- file.path(output_root,paste0(animal,'_components'))
figure_lines <- gsub(paste0('"',animal,'_components"'),encodeString(figure_dir,quote='"'),figure_lines,fixed=TRUE)
# Use a white canvas when saving the supplied PNG plots.
figure_lines <- gsub('dpi = 300', 'dpi = 300, bg = "white"', figure_lines, fixed=TRUE)
eval(parse(text=figure_lines),envir=.GlobalEnv)
expected <- file.path(figure_dir,paste0(animal,'_comp',1:10,'.png'))
stopifnot(all(file.exists(expected)),all(file.info(expected)$size>1000))
if(animal=='A26') {
  eval(parse(text=chunks[['selected protein clusters']]),envir=.GlobalEnv)
  table_path <- file.path(output_root,'Supp_Table_S3_protein_clusters_min10_UPDATED.xlsx')
  sheets <- export_s3(a26.pro.cluster,table_path,min_proteins=10L,expected_clusters=NULL)
  saved_names <- readxl::excel_sheets(table_path)
  stopifnot(identical(saved_names,names(sheets)))
  for(sheet in saved_names) {
    actual <- as.data.frame(readxl::read_excel(table_path,sheet=sheet))
    stopifnot(isTRUE(all.equal(actual,sheets[[sheet]],check.attributes=FALSE,tolerance=1e-10)))
  }
}
dev.off()
writeLines(capture.output(sessionInfo()),file.path(diagnostic_dir,paste0(animal,'_C6_sessionInfo.txt')))
print(summary,row.names=FALSE)
message(animal,': model, plots, and output checks complete at ',format(Sys.time()))
