\name{PeakSegJointHeuristic}
\alias{PeakSegJointHeuristic}
\title{PeakSegJointHeuristic}
\description{Run the PeakSegJoint fast heuristic optimization algorithm, which
gives an approximate solution to a multi-sample Poisson maximum
likelihood segmentation problem. Given S samples, this function
computes a sequence of S+1 PeakSegJoint models, with 0, ..., S
samples with an overlapping peak (maximum of one peak per
sample). This solver runs steps 1-3, and Step3 checks if there are
any more likely models in samples with peak locations which are
the same as all the models detected in Step2. This is guaranteed
as of 24 July 2015 to return a feasible segmentation (seg1 < seg2
> seg3). NB: this function is mostly for internal testing purposes
(search tests/testthat/*.R for 'Heuristic('). For real data use
PeakSegJointSeveral.}
\usage{PeakSegJointHeuristic(profiles, bin.factor = 2L)}
\arguments{
  \item{profiles}{List of data.frames with columns chromStart, chromEnd, count, or
single data.frame with additional column sample.id.}
  \item{bin.factor}{Size of bin pyramid. Bigger values result in slower computation.}
}

\value{List of model fit results, which can be passed to ConvertModelList
for easier interpretation.}

\author{Toby Dylan Hocking}




\examples{

library(PeakSegJoint)
data(H3K36me3.TDH.other.chunk1, envir=environment())
lims <- c(43000000, 43200000) # left
some.counts <-
  subset(H3K36me3.TDH.other.chunk1$counts,
         lims[1] < chromEnd & chromStart < lims[2])
fit <- PeakSegJointHeuristic(some.counts)
converted <- ConvertModelList(fit)
## Normalize profile counts to [0,1].
profile.list <- split(some.counts, some.counts$sample.id)
norm.list <- list()
for(sample.id in names(profile.list)){
  one <- profile.list[[sample.id]]
  max.count <- max(one$count)
  one$count.norm <- one$count/max.count
  norm.list[[sample.id]] <- one
}
norm.df <- do.call(rbind, norm.list)
best.peaks <- transform(converted$peaks, y=peaks*-0.1, what="peaks")

if(interactive() && require(ggplot2)){
  
  ggplot()+
    scale_color_manual(values=c(data="grey50",
                         peaks="deepskyblue",
                         bins="black", segments="green"))+
    geom_step(aes(chromStart/1e3, count.norm, color=what),
              data=data.frame(norm.df, what="data"))+
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

  ggplot(converted$loss, aes(peaks, loss))+
    geom_point()+
    geom_line()

}

}
