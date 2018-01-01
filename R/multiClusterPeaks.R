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
    cluster=as.integer(peaks$chromEnd),
    as.integer(nrow(peaks)),
    PACKAGE="PeakSegJoint")
  peaks$cluster <- res$cluster
  peaks
### peaks data.frame, sorted by chromStart, with an additional column
### cluster.
}, ex=function(){

  library(PeakSegJoint)
  data(chr7.peaks, envir=environment())
  library(ggplot2)
  ggplot()+
    geom_segment(aes(
      chromStart/1e3, sample.id,
      xend=chromEnd/1e3, yend=sample.id),
      data=chr7.peaks)

  clustered <- multiClusterPeaks(chr7.peaks)
  clustered.list <- split(clustered, clustered$cluster)
  clusters.list <- list()
  for(cluster.name in names(clustered.list)){
    clusters.list[[cluster.name]] <- with(
      clustered.list[[cluster.name]], data.frame(
        cluster=cluster[1],
        clusterStart=as.integer(median(chromStart)),
        clusterEnd=as.integer(median(chromEnd))))
  }
  clusters <- do.call(rbind, clusters.list)
  ggplot()+
    geom_segment(aes(
      chromStart/1e3, sample.id,
      color=factor(cluster),
      xend=chromEnd/1e3, yend=sample.id),
      data=clustered)+
  geom_segment(aes(
    clusterStart/1e3, "clusters",
    color=factor(cluster),
    xend=clusterEnd/1e3, yend="clusters"),
    data=clusters)
  
})
