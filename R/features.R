### Compute the feature matrix for this joint segmentation problem.
featureMatrix <- structure(function(profile.list){
  stopifnot(is.list(profile.list))
  stopifnot(is.data.frame(profile.list[[1]]))
  features.by.sample <- list()
  min.chromEnd <- min(sapply(profile.list, with, chromEnd[length(chromEnd)]))
  max.chromStart <- max(sapply(profile.list, with, chromStart[1]))
  bases <- min.chromEnd-max.chromStart
  for(sample.id in names(profile.list)){
    ## Compute feature vector for learning using this segmentation
    ## problem.
    one.list <- profile.list[sample.id]
    sample.counts <- profile.list[[sample.id]]
    loss <- tryCatch({
      fit <- PeakSegJointHeuristic(one.list)
      seg.bases <- with(fit, seg_start_end[2] - seg_start_end[1])
      loss <- sapply(fit$models, "[[", "loss")/seg.bases
      if(!is.finite(loss[2])){
        loss[2] <- loss[1]
      }
      loss
    }, error=function(e){
      bases.vec <- with(sample.counts, chromEnd-chromStart)
      seg.bases <- sum(bases.vec)
      seg.mean <- sum(sample.counts$count  * bases.vec)/seg.bases
      loss0 <- PoissonLoss(sample.counts$count, seg.mean, bases.vec)/seg.bases
      c(loss0, loss0)
    })
    stopifnot(length(loss) == 2)
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
        sum=sum(long),
        loss.0peaks=loss[1],
        loss.1peak=loss[2],
        loss.diff=loss[1]-loss[2])
    suppressWarnings({
      features.by.sample[[sample.id]] <-
        c(feature.vec,
          `log+1`=log(feature.vec+1),
          log=log(feature.vec),
          log.log=log(log(feature.vec)))
    })
  }
  do.call(rbind, features.by.sample)
### Numeric feature matrix (samples x features).
}, ex=function(){
  library(PeakSegJoint)
  data(H3K36me3.TDH.other.chunk1)
  lims <- c(43000000, 43200000) # left
  some.counts <-
    subset(H3K36me3.TDH.other.chunk1$counts,
           lims[1] < chromEnd & chromStart < lims[2])
  profile.list <- ProfileList(some.counts)
  featureMatrix(profile.list)
})
