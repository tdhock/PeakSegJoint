\name{multiClusterPeaks}
\alias{multiClusterPeaks}
\title{multiClusterPeaks}
\description{Cluster peaks into overlapping groups.}
\usage{multiClusterPeaks(peaks)}
\arguments{
  \item{peaks}{data.frame with columns chromStart, chromEnd.}
}

\value{peaks data.frame, sorted by chromStart, with an additional column
cluster.}

\author{Toby Dylan Hocking}




\examples{

library(PeakSegJoint)
data(chr7.peaks)
library(ggplot2)
ggplot()+
  geom_segment(aes(
    chromStart/1e3, sample.id,
    xend=chromEnd/1e3, yend=sample.id),
    data=chr7.peaks)

clustered <- multiClusterPeaks(chr7.peaks)
library(data.table)
clusters <- data.table(clustered)[, list(
  clusterStart=as.integer(median(chromStart)),
  clusterEnd=as.integer(median(chromEnd))
  ), by=cluster]
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

}
