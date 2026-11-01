source('Code/protein_spqc_adjustment.R')
m <- data.frame(pro.sample=c('a','b','c','d','e','f'),
                pro.plate=c('P2','P1','P2','P1','P2','P1'),
                pro.sampletype=c('SPQC','Non-SPQC','Non-SPQC','SPQC','Non-SPQC','SPQC'))
X <- data.frame(p1=c(10,8,14,2,19,4),p2=c(20,21,26,5,31,9))
y <- protein_spqc_adjustment(X,m)
expected <- cbind(p1=c(0,5,4,-1,9,1),p2=c(0,14,6,-2,11,2))
stopifnot(isTRUE(all.equal(unname(y),unname(expected))),
          identical(rownames(y),paste(m$pro.plate,m$pro.sample,sep='::')))
stopifnot(inherits(try(protein_spqc_adjustment(X,m[6:1,]),silent=TRUE),'try-error'))
bad <- m;bad$pro.sampletype[bad$pro.plate=='P2'] <- 'Non-SPQC'
stopifnot(inherits(try(protein_spqc_adjustment(X,bad),silent=TRUE),'try-error'))
# Repeated SPQC labels are valid and must retain distinct row identities.
dup_m <- rbind(m,m[1,]);dup_x <- rbind(X,X[1,])
rownames(dup_m) <- rownames(dup_x) <- NULL
dup_y <- protein_spqc_adjustment(dup_x,dup_m)
stopifnot(nrow(dup_y)==7L,!anyDuplicated(rownames(dup_y)),all(dup_y[7,]==0))
cat('Interleaved-plate SPQC adjustment, sample alignment, and missing-QC checks passed.\n')
