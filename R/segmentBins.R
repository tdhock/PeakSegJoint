segmentBins <- function
### Run binSum and PeakSegDP on a chromosome subset.
(compressed,
### data.frame with columns chromStart, chromEnd, count of coverage on
### one chromosome.
 bin.chromStart,
### integer: chromStart of first bin.
 bin.size,
### integer: bases per bin.
 n.bins=2000L,
### integer: number of bins.
 maxPeaks=9L
### integer: maximum number of peaks in segmentation model.
 ){
  binSum.seconds <- system.time({
    bins <- 
      binSum(compressed,
             bin.chromStart=bin.chromStart,
             bin.size=bin.size,
             n.bins=n.bins)
  })[["elapsed"]]

  maxSegments <- nrow(bins)
  this.max <- as.integer(min((maxSegments-1)/2, maxPeaks))
  cDPA.seconds <- system.time(suppressWarnings({
    fit <- PeakSegDP(bins, maxPeaks = this.max)
  }))[["elapsed"]]

  ## Compute feature vector for learning using this segmentation
  ## problem.
  bases <- with(bins, chromEnd-chromStart)
  long <- rep(bins$count, bases)
  n.bases <- sum(bases)
  n.data <- nrow(bins)

  uq <- quantile(bins$count)
  chrom.quartile <- quantile(compressed$count)
  feature.vec <-
    c(unweighted.quartile=uq,
      weighted.quartile=quantile(long),
      chrom.quantile.ratio=uq/chrom.quartile,
      chrom.quantile=chrom.quartile,
      unweighted.mean=mean(bins$count),
      weighted.mean=mean(long),
      bases=n.bases,
      data=n.data)
  
  suppressWarnings({
    features <-
      c(feature.vec,
        `log+1`=log(feature.vec+1),
        log=log(feature.vec),
        log.log=log(log(feature.vec)))
  })

  all.loss <- data.frame(fit$error)
  all.loss$cummin <- cummin(all.loss$error)
  loss <- subset(all.loss, error == cummin)
  rownames(loss) <- loss$segments
  bases <- with(fit$segments[1,], chromEnd-chromStart)
  in.sqrt <- 1.1 + log(bases / loss$segments)
  in.square <- 1 + 4 * sqrt(in.sqrt)
  complexity <- in.square * in.square * loss$segments
  
  exact <- with(loss, exactModelSelection(error, complexity, peaks))

  list(features=features, bins=bins, modelSelection=exact, fit=fit,
       binSum.seconds=binSum.seconds, cDPA.seconds=cDPA.seconds)
### List of features, bins, modelSelection, fit, and timings.
}
