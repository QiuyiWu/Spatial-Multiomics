# Rscript Code/regenerate_QCFix.R FRESH_RUN_DIRECTORY FIGURE_OUTPUT_DIRECTORY
# Regenerates all three QCFix corrections and the complete manuscript figure set.
# Use a restored R environment. PROTEINFIX_RSCRIPT_FLAGS=--vanilla is an explicit
# alternative for an independently configured library (record sessionInfo).
args <- commandArgs(TRUE);stopifnot(length(args)==2L)
here::i_am('README.md')
run <- normalizePath(args[1],mustWork=FALSE)
out <- normalizePath(args[2],mustWork=FALSE)
if(dir.exists(run) && length(list.files(run,all.files=TRUE,no..=TRUE)))
  stop('Use a fresh empty run directory to prevent stale intermediate reuse')
if(dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE)))
  stop('Use a fresh empty figure directory to prevent mixing old and new exports')
dir.create(run,recursive=TRUE,showWarnings=FALSE)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
logs <- file.path(run,'logs');dir.create(logs)
Sys.setenv(PROTEINFIX_MODE='reproduce',PROTEINFIX_RUN_DIR=run,
           PROTEINFIX_C6_DIAGNOSTIC_DIR=file.path(logs,'C6'))
flags <- Sys.getenv('PROTEINFIX_RSCRIPT_FLAGS','')
flags <- if(nzchar(flags))strsplit(flags,' ',fixed=TRUE)[[1]] else character()
execute <- function(script,arguments=character(),label=basename(script)) {
 log <- file.path(logs,paste0(label,'.log'))
 status <- system2(file.path(R.home('bin'),'Rscript'),
                   c(flags,shQuote(script),vapply(arguments,shQuote,character(1))),
                   stdout=log,stderr=log)
 if(status!=0L)stop('Stage failed: ',label,'; inspect ',log)
 message('Completed ',label)
}
# Notebook and selection stages run in separate sessions to avoid stale objects.
stage <- file.path(run,'preprocess_and_select.R')
writeLines(c("here::i_am('README.md')",
 "run <- Sys.getenv('PROTEINFIX_RUN_DIR')",
 "pdf(file.path(run,'logs','A27_preprocessing.pdf'))",
 "script <- file.path(run,'logs','A27_met_executed.R')",
 "knitr::purl('Code/A27/Analysis_BatchEffect_A27_met_adj.Rmd',output=script,documentation=0,quiet=TRUE)",
 "source(script)",
 "write.csv(data.frame(level=names(plasma_qc_counts),removed=as.integer(plasma_qc_counts)),file.path(run,'plasma_qc_exclusions.csv'),row.names=FALSE)",
 "source('Code/A27/HeldoutPreprocessing.R')",
 "source('Code/A27/CombineClustering.R')",
 "dev.off()"),stage)
execute(stage,label='A27_preprocess_select_cluster')
for(animal in c('A26','A27'))execute('Code/regenerate_C6_outputs.R',c(animal,out),paste0(animal,'_C6'))
execute('Code/regenerate_projection.R',file.path(out,'A27_to_A26_components'),'projection')
pca <- file.path(run,'pca.R')
writeLines(c(sprintf('options(proteinfix.pca_export_dir=%s)',encodeString(out,quote='"')),
             "source('Code/PCA_Figure.R')"),pca)
execute(pca,label='PCA')
execute('Code/regenerate_Figure4.R',file.path(out,'Figure4'),'Figure4')
execute('Code/regenerate_FigureS5.R',file.path(out,'FigureS5'),'FigureS5')
execute('Code/regenerate_Figure5.R',file.path(out,'Figure5'),'Figure5')
execute('Code/assemble_manuscript_components.R',c(out,out),'component_assembly')
exports <- c('post_PCA_4panel_shared_legend.tiff'='Figure2.tiff',
             'pre_PCA_4panel_shared_legend.tiff'='FigureS2.tiff',
             'Figure4/Figure4_corrected.tiff'='Figure4.tiff',
             'Figure5/fig5_pca.tiff'='Figure5.tiff',
             'FigureS5/FigureS5_corrected.tiff'='FigureS5.tiff')
for(source in names(exports))stopifnot(file.copy(file.path(out,source),file.path(out,exports[[source]]),overwrite=TRUE))
execute('Code/validate_QCFix.R',run,'validation')
execute('Code/compare_reproduction.R',run,'comparison')
writeLines(capture.output(sessionInfo()),file.path(run,'sessionInfo.txt'))
message('Finished. Review regenerated Data/ and comparison.csv before promoting references.')
