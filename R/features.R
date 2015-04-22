### Compute the feature matrix for this joint segmentation problem.
featureMatrix <- function(profile.list){
  features.by.sample <- list()
  min.chromEnd <- min(sapply(profile.list, with, chromEnd[length(chromEnd)]))
  max.chromStart <- max(sapply(profile.list, with, chromStart[1]))
  bases <- min.chromEnd-max.chromStart
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
  do.call(rbind, features.by.sample)
}
