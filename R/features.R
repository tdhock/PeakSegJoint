### Compute the feature matrix for this joint segmentation problem.
computeFeatures <- function(){
  features.by.sample <- list()
  for(sample.id in names(profile.list)){
    ## Compute feature vector for learning using this segmentation
    ## problem.
    sample.counts <- profile.list[[sample.id]]
    too.long <- with(sample.counts, rep(count, chromEnd-chromStart))
    too.long.pos <- with(sample.counts, {
      (chromStart[1]+1):chromEnd[length(chromEnd)]
    })
    stopifnot(length(too.long) == length(too.long.pos))
    keep <- max.chromStart < too.long.pos & too.long.pos <= min.chromEnd
    long <- as.numeric(too.long[keep])
    stopifnot(length(long) == bases)
    feature.vec <-
      c(quartile=quantile(long),
        mean=mean(long),
        sd=sd(long),
        mad=mad(long),
        bases=bases,
        sum=sum(long))
    suppressWarnings({
      features.by.sample[[sample.id]] <-
        c(feature.vec,
          `log+1`=log(feature.vec+1),
          log=log(feature.vec),
          log.log=log(log(feature.vec)))
    })
  }
  feature.mat <- do.call(rbind, features.by.sample)
}
