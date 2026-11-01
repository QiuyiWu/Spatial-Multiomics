# Rscript Code/regenerate_projection.R OUTPUT_COMPONENT_DIRECTORY
args <- commandArgs(TRUE);stopifnot(length(args)==1L)
here::i_am('README.md')
out <- normalizePath(args[1],mustWork=FALSE)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
options(proteinfix.projection_export_dir=out)
script <- file.path(out,'projection_executed.R')
knitr::purl('Code/A26A27_analysis_sCCA.Rmd',output=script,documentation=0,quiet=TRUE)
pdf(file.path(out,'projection_diagnostics.pdf'))
source(script,echo=FALSE)
dev.off()
saved <- readRDS(proteinfix_input('A27_to_A26_projection.rds'))
stopifnot(isTRUE(all.equal(saved$protein_scores,X26_pro_common %*% w27_pro_common,check.attributes=FALSE,tolerance=1e-12)),
          isTRUE(all.equal(saved$metabolite_scores,X26_met_common %*% w27_met_common,check.attributes=FALSE,tolerance=1e-12)))
writeLines(c(paste(length(common_pro),'shared proteins;',length(common_met),'shared metabolites.'),
             'All component PNGs and individual TIFFs use the stored fitted-model signs.',
             'No presentation-only sign flips are applied.'),file.path(out,'projection_orientation.txt'))
cat('RESULT shared proteins',length(common_pro),'metabolites',length(common_met),'\n')
