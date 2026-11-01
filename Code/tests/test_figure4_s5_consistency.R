# Rscript Code/tests/test_figure4_s5_consistency.R FIGURE4_DIRECTORY S5_DIRECTORY [PREVIOUS_FIGURE4_DIRECTORY]
args <- commandArgs(TRUE);stopifnot(length(args)%in%c(2L,3L))
for(animal in c('A26','A27')) {
  f <- readRDS(file.path(args[1],paste0(animal,'_Figure4_results.rds')))
  s <- readRDS(file.path(args[2],paste0(animal,'_FigureS5_results.rds')))
  stopifnot(f$protein_representation=='loading_weighted_clusters',
            f$protein_clusters==s$clusters,f$n==s$n,identical(f$orders,s$orders),
            isTRUE(all.equal(f$pair_correlations$met,s$met,tolerance=1e-12)),
            isTRUE(all.equal(f$pair_correlations$pro,s$procluster,tolerance=1e-12)))
  if(length(args)==3L) {
    old <- readRDS(file.path(args[3],paste0(animal,'_Figure4_results.rds')))
    stopifnot(isTRUE(all.equal(f$correlations$met,old$correlations$met,tolerance=1e-12)))
  }
  cat(animal,': Figure 4/S5 pair correlations match;',f$protein_clusters,'protein clusters.\n')
}
