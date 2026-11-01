# Shared routing: reference inputs remain in Data/; regenerated artifacts are isolated.
proteinfix_mode <- function() {
  mode <- getOption('proteinfix.mode', Sys.getenv('PROTEINFIX_MODE', 'reproduce'))
  if (!mode %in% c('reproduce','update')) stop('proteinfix.mode must be reproduce or update')
  mode
}
proteinfix_run_dir <- function() {
  path <- getOption('proteinfix.run_dir', Sys.getenv('PROTEINFIX_RUN_DIR', here::here('outputs','reproduction')))
  dir.create(path,recursive=TRUE,showWarnings=FALSE)
  path <- normalizePath(path,mustWork=TRUE)
  reference <- normalizePath(here::here('Data'),mustWork=TRUE)
  if (identical(path,reference) || startsWith(paste0(path,'/'),paste0(reference,'/')))
    stop('The reproduction directory must be outside Data/')
  candidate <- file.path(path,'Data')
  if(dir.exists(candidate) && identical(normalizePath(candidate),reference))
    stop('The reproduction Data directory must not be the reference Data directory')
  path
}
proteinfix_output <- function(name) {
  stopifnot(length(name)==1L, basename(name)==name)
  folder <- if(proteinfix_mode()=='update') here::here('Data') else file.path(proteinfix_run_dir(),'Data')
  dir.create(folder,recursive=TRUE,showWarnings=FALSE)
  path <- file.path(folder,name)
  message('Writing ',proteinfix_mode(),' output: ',path)
  path
}
proteinfix_input <- function(name) {
  stopifnot(length(name)==1L, basename(name)==name)
  path <- here::here('Data',name)
  if(proteinfix_mode()=='reproduce') {
    candidate <- file.path(proteinfix_run_dir(),'Data',name)
    if(file.exists(candidate)) path <- candidate
  }
  if(!file.exists(path)) stop('Missing input: ',path)
  message('Reading input: ',path)
  path
}
proteinfix_graph <- function(name, scores, xy) {
  if(!identical(rownames(scores),rownames(xy)) || !all(is.finite(as.matrix(xy))) || anyDuplicated(rownames(xy)))
    stop('Spatial graph inputs must have aligned, unique sample IDs and finite coordinates')
  rebuild <- isTRUE(getOption('proteinfix.rebuild_graph',identical(tolower(Sys.getenv('PROTEINFIX_REBUILD_GRAPH','false')),'true')))
  if(rebuild) {
    graph <- adegenet::spca(x=scores,xy=xy,center=TRUE,scale=TRUE,scannf=FALSE,
                           nfposi=3,nfnega=0,type=1,ask=FALSE,plot.nb=FALSE)
  } else {
    e <- new.env(parent=emptyenv())
    load(proteinfix_input(name),envir=e)
    graph <- e$spca_delaunay_metCluster
  }
  if(is.null(graph) || !identical(rownames(graph$li),rownames(scores)) ||
     !identical(rownames(graph$xy),rownames(xy)) ||
     !isTRUE(all.equal(unname(as.matrix(graph$xy)),unname(as.matrix(xy)),tolerance=1e-12)))
    stop('Saved graph does not match the current sample order/coordinates. Explicitly set options(proteinfix.rebuild_graph=TRUE) to rebuild.')
  n <- nrow(scores)
  A <- matrix(0,n,n)
  if(length(graph$lw$neighbours)!=n) stop('Saved graph has the wrong number of nodes')
  for(i in seq_len(n)) {
    j <- graph$lw$neighbours[[i]]
    if(any(j<1L | j>n | j==i)) stop('Invalid saved graph neighbor index')
    A[i,j] <- 1
  }
  if(!isSymmetric(A) || any(rowSums(A)==0)) stop('Saved graph must be symmetric and have no isolated samples')
  if(rebuild) {
    spca_delaunay_metCluster <- graph
    save(spca_delaunay_metCluster,file=proteinfix_output(name))
  }
  graph
}
