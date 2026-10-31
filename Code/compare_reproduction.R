# Compare regenerated analytical contents, not serialization/Excel timestamps.
# Run from the repository root: Rscript Code/compare_reproduction.R [run_directory]
source(here::here('Code','pipeline_paths.R'))
args <- commandArgs(TRUE)
if(length(args)) options(proteinfix.run_dir=args[1])
root <- proteinfix_run_dir()
files <- list.files(file.path(root,'Data'),full.names=TRUE)
compare_file <- function(path) {
  reference <- here::here('Data',basename(path))
  if(!file.exists(reference)) return('No reference file')
  ext <- tolower(tools::file_ext(path))
  read_data <- function(p) {
    if(ext=='xlsx') {
      sheets <- readxl::excel_sheets(p)
      return(setNames(lapply(sheets,function(s)as.data.frame(readxl::read_xlsx(p,sheet=s))),sheets))
    }
    if(ext=='rdata') {
      e <- new.env(parent=emptyenv());load(p,e)
      objects <- mget(sort(ls(e)),e)
      if('spca_delaunay_metCluster' %in% names(objects)) {
        g <- objects$spca_delaunay_metCluster
        return(list(samples=rownames(g$li),xy=g$xy,neighbors=g$lw$neighbours))
      }
      return(objects)
    }
    if(ext=='rds') return(readRDS(p))
    stop('Unsupported comparison format: ',ext)
  }
  answer <- all.equal(read_data(reference),read_data(path),tolerance=1e-8,check.attributes=TRUE)
  if(isTRUE(answer)) 'Match' else paste(answer,collapse='; ')
}
if(!length(files)) stop('No regenerated Data files found in ',root)
report <- data.frame(file=basename(files),result=vapply(files,function(p)
  tryCatch(compare_file(p),error=function(e)paste('Comparison failed:',conditionMessage(e))),character(1)))
write.csv(report,file.path(root,'comparison.csv'),row.names=FALSE)
print(report,row.names=FALSE)
if(any(report$result!='Match')) warning('Some regenerated outputs differ; inspect comparison.csv before updating references.')
