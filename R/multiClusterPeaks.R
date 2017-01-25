multiClusterPeaks <- structure(function
### Cluster peaks into overlapping groups.
(peaks
### data.frame with columns chromStart, chromEnd.
 ){
  stopifnot(is.data.frame(peaks))
  stopifnot(is.numeric(peaks$chromStart))
  stopifnot(is.numeric(peaks$chromEnd))
  stopifnot(nrow(peaks) > 0)
  peaks <- data.frame(peaks)[order(peaks$chromStart), ]
  with(peaks, stopifnot(chromStart < chromEnd))
  sample.vec <- if(is.character(peaks$sample.id)){
    factor(peaks$sample.id)
  }else if(is.factor(peaks$sample.id) || is.integer(peaks$sample.id)){
    peaks$sample.id
  }else stop("sample.id column must be character, factor, or integer")
  res <- .C(
    "multiClusterPeaks_interface",
    as.integer(peaks$chromStart),
    as.integer(peaks$chromEnd),
    as.integer(sample.vec),
    as.integer(nrow(peaks)),
    cluster=as.integer(peaks$chromEnd),
    PACKAGE="PeakSegJoint")
  peaks$cluster <- res$cluster
  peaks
### peaks data.frame, sorted by chromStart, with an additional column
### cluster.
}, ex=function(){

  library(PeakSegJoint)
  data(chr7.peaks)
  library(ggplot2)
  ggplot()+
    geom_segment(aes(
      chromStart/1e3, sample.id,
      xend=chromEnd/1e3, yend=sample.id),
      data=chr7.peaks)

  clustered <- multiClusterPeaks(chr7.peaks)
  ggplot()+
    geom_segment(aes(
      chromStart/1e3, sample.id,
      color=factor(cluster),
      xend=chromEnd/1e3, yend=sample.id),
      data=clustered)
  
  gg+geom_segment(aes(
    clusterStart/1e3, "clusters",
    xend=clusterEnd/1e3, yend="clusters"), data=clusters)
  
})
