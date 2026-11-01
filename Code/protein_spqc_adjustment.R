# Stable row identities include occurrence suffixes for legitimate repeated SPQC labels.
protein_sample_keys <- function(metadata) {
  plate <- as.character(metadata$pro.plate)
  sample <- as.character(metadata$pro.sample)
  type <- as.character(metadata$pro.sampletype)
  if(anyNA(plate) || anyNA(sample) || anyNA(type) ||
     any(!nzchar(plate)) || any(!nzchar(sample)))
    stop('SPQC adjustment requires complete sample, plate, and sample-type metadata')
  make.unique(paste(plate,sample,sep='::'))
}

# Subtract plate-specific SPQC means without changing sample order.
protein_spqc_adjustment <- function(values, metadata) {
  X <- as.matrix(values)
  stopifnot(is.numeric(X), nrow(X)==nrow(metadata), all(is.finite(X)),
            identical(rownames(values),rownames(metadata)))
  plate <- as.character(metadata$pro.plate)
  sample <- as.character(metadata$pro.sample)
  type <- as.character(metadata$pro.sampletype)
  keys <- protein_sample_keys(metadata)
  rownames(X) <- keys
  adjusted <- X
  for(p in unique(plate)) {
    idx <- which(plate==p)
    qc <- idx[type[idx]=='SPQC']
    if(!length(qc)) stop('No SPQC rows for plate ',p)
    adjusted[idx,] <- sweep(X[idx,,drop=FALSE],2,
                           colMeans(X[qc,,drop=FALSE]),'-')
  }
  stopifnot(identical(rownames(adjusted),keys),all(is.finite(adjusted)))
  adjusted
}
