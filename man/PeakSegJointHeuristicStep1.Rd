\name{PeakSegJointHeuristicStep1}
\alias{PeakSegJointHeuristicStep1}
\title{PeakSegJointHeuristicStep1}
\description{Run the first step of the PeakSegJoint fast heuristic optimization
algorithm. This is the GridSearch subroutine of the JointZoom algorithm
in arXiv:1506.01286. NB: this function is only for testing the C code
against the R implementation (search tests/testthat/*.R for Step1).
For real data see
\code{\link{PeakSegJointSeveral}}.}
\usage{PeakSegJointHeuristicStep1(profiles, 
    bin.factor = 2L)}
\arguments{
  \item{profiles}{List of data.frames with columns chromStart, chromEnd, count, or
single data.frame with additional column sample.id.}
  \item{bin.factor}{Size of bin pyramid. Bigger values result in slower computation.}
}

\value{List of model fit results, which can be passed to \code{\link{ConvertModelList}}
for easier interpretation.}

\author{Toby Dylan Hocking <toby.hocking@r-project.org> [aut, cre]}




\examples{
library(PeakSegJoint)
library(ggplot2)
data(H3K36me3.TDH.other.chunk1, envir=environment())
lims <- c(43000000, 43200000) # left
some.counts <-
  subset(H3K36me3.TDH.other.chunk1$counts,
         lims[1] < chromEnd & chromStart < lims[2])
fit <- PeakSegJointHeuristicStep1(some.counts)
## Compute bins to show on plot as well.
profile.list <- split(some.counts, some.counts$sample.id)
bin.list <- list()
norm.list <- list()
for(sample.id in names(profile.list)){
  one <- profile.list[[sample.id]]
  max.count <- max(one$count)
  bins <- binSum(one, fit$bin_start_end[1], fit$bases_per_bin, fit$n_bins)
  stopifnot(fit$n_bins == nrow(bins))
  bins$mean <- with(bins, count/(chromEnd-chromStart))
  bins$mean.norm <- bins$mean/max.count
  stopifnot(bins$count >= 0)
  one$count.norm <- one$count/max.count
  norm.list[[sample.id]] <- one
  bin.list[[sample.id]] <- data.frame(sample.id, bins)
}
bin.df <- do.call(rbind, bin.list)
norm.df <- do.call(rbind, norm.list)
converted <- ConvertModelList(fit)
best.peaks <- transform(converted$peaks, y=peaks*-0.1, what="peaks")

if(require(ggplot2) && interactive()){

  ggplot()+
    scale_color_manual(values=c(data="grey50",
                                peaks="deepskyblue",
                                bins="black", segments="green"))+
    geom_step(aes(chromStart/1e3, count.norm, color=what),
              data=data.frame(norm.df, what="data"))+
    geom_segment(aes(chromStart/1e3, mean.norm,
                     xend=chromEnd/1e3, yend=mean.norm,
                     color=what),
                 data=data.frame(bin.df, what="bins"))+
    geom_segment(aes(chromStart/1e3, y,
                     xend=chromEnd/1e3, yend=y,
                     color=what),
                 size=1,
                 data=best.peaks)+
    geom_text(aes(chromStart/1e3, y,
                  label=paste0(peaks, " peak",
                               ifelse(peaks==1, "", "s"), " "),
                  color=what),
              hjust=1,
              size=3,
              vjust=0.5,
              data=best.peaks)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., scales="free")

}

}
