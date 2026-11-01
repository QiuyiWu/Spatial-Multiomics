# Rscript Code/validate_QCFix.R RUN_DIRECTORY
args <- commandArgs(TRUE);stopifnot(length(args)==1L)
here::i_am('README.md')
options(proteinfix.mode='reproduce',proteinfix.run_dir=args[1])
source('Code/pipeline_paths.R')
load(proteinfix_input('met_a27_PCAdata.RData'))
stopifnot(identical(dim(pre_met_a27),c(581L,300L)),
          identical(dim(post_met_a27),c(548L,300L)),
          !any(grepl('QC Level',pre_met_a27_metadata$sample_id)),
          sum(pre_met_a27_metadata$sample_type=='SPQC')==33L,
          all(is.finite(pre_met_a27)),all(is.finite(post_met_a27)))
load(proteinfix_input('a27.RData'))
load(proteinfix_input('A27_PChclust.RData'))
stopifnot(nrow(a27.met.data)==527L,ncol(a27.met.data)-2L==300L,
          identical(rownames(a27.met.cluster$score),rownames(a27.pro.cluster$score)),
          setequal(rownames(a27.met.cluster$loading),make.names(colnames(a27.met.data)[-c(1,2)])),
          all(is.finite(a27.met.cluster$loading)),all(is.finite(a27.pro.cluster$loading)))
old <- new.env();load('Data/A27_PChclust.RData',old)
cat('Protein clustering unchanged:',isTRUE(all.equal(old$a27.pro.cluster,a27.pro.cluster,tolerance=1e-10)),'\n')
load(proteinfix_input('a27_score.RData'))
stopifnot(identical(dim(scca_res$u),c(ncol(a27.pro.cluster$score),10L)),
          identical(dim(scca_res$v),c(ncol(a27.met.cluster$score),10L)),
          all(is.finite(scca_res$u)),all(is.finite(scca_res$v)))
projection <- readRDS(proteinfix_input('A27_to_A26_projection.rds'))
stopifnot(length(projection$shared_proteins)==9250L, length(projection$shared_metabolites)==277L)
stopifnot(identical(dim(projection$protein_scores),c(267L,10L)),
          identical(dim(projection$metabolite_scores),c(267L,10L)),
          all(is.finite(projection$protein_scores)),all(is.finite(projection$metabolite_scores)))
f5 <- readRDS(proteinfix_input('Figure5_PCA_results.rds'))
stopifnot(nrow(f5$clean_data)==1013L,nrow(f5$pca$rotation)==252L,
          sum(f5$outlier_audit$flagged)==14L,
          sum(f5$outlier_audit$legacy_removed)==21L,
          sum(f5$outlier_audit$extra_legacy_removal)==7L)
for(m in f5$matches)stopifnot(identical(m$pc1$Feature,rownames(m$weights)),
                            !anyDuplicated(m$mapping$key),all(is.finite(m$weights)))
cat('All QCFix output checks passed.\n')
