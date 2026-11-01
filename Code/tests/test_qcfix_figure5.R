source('Code/figure5_checks.R')
d <- data.frame(Subject=c('A26','A27','U5'),Sample=c('123','123','9'))
s <- transform(d,outlier=c(TRUE,FALSE,FALSE))
a <- figure5_outliers(d,s)
stopifnot(sum(a$flagged)==1L,sum(a$legacy_removed)==2L,
          identical(which(a$extra_legacy_removal),2L))
stopifnot(inherits(try(figure5_outliers(d,s[3:1,]),silent=TRUE),'try-error'))
p <- data.frame(Feature=c('PC aa C36:4','Ala','not shared'),PC1_Loading=c(.2,.3,.4))
w <- matrix(1:6,3,dimnames=list(make.names(c('PC aa C36:4','Ala','not shared')),c('Comp1','Comp2')))
m <- figure5_match_loadings(p,w[3:1,,drop=FALSE])
stopifnot(nrow(m$pc1)==3L,m$literal_matches==1L,
          identical(m$pc1$Feature,rownames(m$weights)),
          identical(m$pc1$OriginalFeature,p$Feature),
          identical(unname(m$weights[,1]),1:3))
p$Feature <- c('same name','same.name','other')
stopifnot(inherits(try(figure5_match_loadings(p,w),silent=TRUE),'try-error'))
cat('Figure 5 repeated-ID removal, alignment, name normalization, and collision checks passed.\n')
