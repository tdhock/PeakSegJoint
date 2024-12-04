### Compute the feature matrix for this joint segmentation problem.
featureMatrixJoint <- structure(function(profile.list){
  chromEnd <- chromStart <- NULL
  stopifnot(is.list(profile.list))
  stopifnot(is.data.frame(profile.list[[1]]))
  features.by.sample <- list()
  ## We would get the following when one sample has no coverage rows:
  ## Warning messages:
  ## 1: In max(chromEnd) : no non-missing arguments to max; returning -Inf
  ## 2: In min(chromStart) : no non-missing arguments to min; returning Inf
  suppressWarnings({
    sample.chromEnds <- sapply(profile.list, with, max(chromEnd))
    sample.chromStarts <- sapply(profile.list, with, min(chromStart))
  })
  min.chromEnd <- as.integer(max(sample.chromEnds))
  max.chromStart <- as.integer(min(sample.chromStarts))
  bases <- min.chromEnd-max.chromStart
  for(sample.id in names(profile.list)){
    ## Compute feature vector for learning using this segmentation
    ## problem.
    one.list <- profile.list[sample.id]
    sample.counts <- profile.list[[sample.id]]
    loss <- tryCatch({
      fit <- PeakSegJointHeuristic(one.list)
      seg.bases <- with(fit, bin_start_end[2] - bin_start_end[1])
      loss <- sapply(fit$models, "[[", "loss")/seg.bases
      if(!is.finite(loss[2])){
        loss[2] <- loss[1]
      }
      loss
    }, error=function(e){
      if(!identical(e$message, "bin factor too large")){
        stop(e)
      }else{
        bases.vec <- with(sample.counts, chromEnd-chromStart)
        seg.bases <- sum(bases.vec)
        seg.mean <- sum(sample.counts$count  * bases.vec)/seg.bases
        loss0 <- PoissonLoss(sample.counts$count, seg.mean, bases.vec)/seg.bases
        c(loss0, loss0)
      }
    })
    stopifnot(length(loss) == 2)
    bins <-
      binSum(sample.counts, max.chromStart,
             bin.size=1L, n.bins=bases,
             empty.as.zero=TRUE)
    long <- bins$count
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
  data(H3K36me3.TDH.other.chunk1, envir=environment())
  lims <- c(43000000, 43200000) # left
  some.counts <-
    subset(H3K36me3.TDH.other.chunk1$counts,
           lims[1] < chromEnd & chromStart < lims[2])
  profile.list <- ProfileList(some.counts)
  featureMatrixJoint(profile.list)
})
