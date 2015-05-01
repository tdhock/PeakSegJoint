clusterPeaks <- structure(function
### Cluster peaks into overlapping groups.
(peaks
### data.frame with columns chromStart, chromEnd.
 ){
  stopifnot(is.data.frame(peaks))
  stopifnot(is.numeric(peaks$chromStart))
  stopifnot(is.numeric(peaks$chromEnd))
  stopifnot(nrow(peaks) > 0)
  peaks <- peaks[order(peaks$chromStart), ]
  res <- .C("clusterPeaks_interface",
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
  unordered <-
    data.frame(chromStart=c(11, 12, 1, 2, 3, 6, 7),
               chromEnd=c(13, 14, 5, 8, 4, 9, 10),
               sample.id=factor(paste0("sample.", c(1, 2, 1, 2, 3, 4, 5))))
  clustered <- clusterPeaks(unordered)
  library(ggplot2)
  bases <- geom_text(aes(position, "base", label=base),
                     data=data.frame(position=1:20, base="N"))
  gg <- ggplot()+bases+ylab("")+
    scale_x_continuous(breaks=1:20)
  
  gg+
    geom_segment(aes(chromStart+1/2, sample.id,
                     xend=chromEnd+1/2, yend=sample.id),
                     data=clustered)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(.~cluster, labeller=label_both, scales="free", space="free")

  gg+
    geom_segment(aes(chromStart+1/2, sample.id,
                     xend=chromEnd+1/2, yend=sample.id, color=factor(cluster)),
                     data=clustered)
})
