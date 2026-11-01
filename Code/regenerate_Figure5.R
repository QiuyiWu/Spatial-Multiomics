# Rscript Code/regenerate_Figure5.R OUTPUT_DIRECTORY
# Runs the notebook verbatim, including its PI-confirmed QC3 preprocessing policy.
args <- commandArgs(TRUE)
stopifnot(length(args)==1L)
here::i_am('README.md')
out <- normalizePath(args[1],mustWork=FALSE)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
options(proteinfix.figure5_export_dir=out)
script <- file.path(out,'Figure5_executed.R')
knitr::purl('Code/loadingComp_U5A26A27_met.Rmd',output=script,documentation=0,quiet=TRUE)
pdf(file.path(out,'Figure5_diagnostics.pdf'))
source(script,echo=FALSE)
dev.off()
stopifnot(all(is.finite(met_combined_mat)),
          nrow(met_data.comb.clean)==sum(!figure5_outlier_audit$flagged))
write.csv(figure5_outlier_audit,file.path(out,'Figure5_outlier_audit.csv'),row.names=FALSE)
write.csv(pca_data.comb,file.path(out,'Figure5_PCA_scores.csv'),row.names=FALSE)
write.csv(pc1_loadings_df,file.path(out,'Figure5_PC1_loadings.csv'),row.names=FALSE)
write.csv(cor_results.a26,file.path(out,'Figure5_A26_loading_correlations.csv'),row.names=FALSE)
write.csv(cor_results.a27,file.path(out,'Figure5_A27_loading_correlations.csv'),row.names=FALSE)
for(animal in c('a26','a27')) {
  m <- get(paste0('matched_',animal)); res <- get(paste0('res_',animal))
  write.csv(m$mapping,file.path(out,paste0('Figure5_',animal,'_feature_matches.csv')),row.names=FALSE)
  writeLines(m$unmatched,file.path(out,paste0('Figure5_',animal,'_unmatched.txt')))
  write.csv(res$overlap_df,file.path(out,paste0('Figure5_',animal,'_overlap.csv')),row.names=FALSE)
  ggplot2::ggsave(file.path(out,paste0('Figure5_',animal,'_overlap.png')),res$plot,width=9,height=6,dpi=180)
}
saveRDS(list(pca=pca_combined,clean_data=met_data.comb.clean,
             outlier_audit=figure5_outlier_audit,
             matches=list(A26=matched_a26,A27=matched_a27),
             correlations=list(A26=cor_results.a26,A27=cor_results.a27)),
        file.path(out,'Figure5_PCA_results.rds'))
file.copy(file.path(out,'Figure5_PCA_results.rds'), proteinfix_output('Figure5_PCA_results.rds'), overwrite=TRUE)
ggplot2::ggsave(file.path(out,'Figure5_preview.png'),p_fig5,width=fig_w_mm,height=fig_h_mm,units='mm',dpi=180,bg='white')
writeLines(capture.output(sessionInfo()),file.path(out,'Figure5_sessionInfo.txt'))
cat('RESULT: features',ncol(met_combined_mat),'samples',nrow(met_combined_mat),
    'flagged',sum(figure5_outlier_audit$flagged),'legacy removed',sum(figure5_outlier_audit$legacy_removed),'\n')
print(table(met_data.comb.clean$Subject))
cat('PC1/PC2',pc1_var,pc2_var,'matches',nrow(matched_a26$pc1),nrow(matched_a27$pc1),'\n')
