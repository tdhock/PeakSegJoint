\name{PeakSegJointSeveral}
\alias{PeakSegJointSeveral}
\title{PeakSegJointSeveral}
\description{Run the PeakSegJoint heuristic segmentation algorithm with several
different bin.factor values, keeping only the models with lowest
Poisson loss for each peak size. This algorithm gives an
approximate solution to the following multi-sample constrained
maximum likelihood segmentation problem. If there are S samples
total, we look for the most likely common peak in $s\in{0, ..., S}$
samples. We solve the equivalent minimization problem using the
Poisson loss seg.mean - count.data * log(seg.mean), from the first
base to the last base of profiles. The optimization variables are
the segment means, of which there can be either 1 value (no peak)
or 3 values (peak) in each sample. If there are 3 segments then
two constraints are applied: (1) the changes in mean must occur at
the same position in each sample, and (2) the changes must be up
and then down (mean1 < mean2 > mean3).}
\usage{PeakSegJointSeveral(profiles, bin.factors = 2:7)}
\arguments{
  \item{profiles}{data.frame or list of them from ProfileList.}
  \item{bin.factors}{integer vector of optimization parameters >= 2. Larger values are
slower. Using more values is slower since it tells the algorithm
to search more of the model space, yielding solution which is
closer to the global optimum.}
}

\value{List of model fit results, which can be passed to ConvertModelList
for easier interpretation.}

\author{Toby Dylan Hocking}




\examples{
data(H3K4me3.TDH.other.chunk8)
bf.vec <- c(2, 3, 5)
fit.list <-
  list(several=PeakSegJointSeveral(H3K4me3.TDH.other.chunk8, bf.vec))
for(bf in bf.vec){
  fit.list[[paste(bf)]] <-
    PeakSegJointHeuristicStep2(H3K4me3.TDH.other.chunk8, bf)
}
loss.list <- list()
segs.by.peaks.fit <- list()
for(fit.name in names(fit.list)){
  fit <- fit.list[[fit.name]]
  loss.list[[fit.name]] <- sapply(fit$models, "[[", "loss")
  converted <- ConvertModelList(fit)
  segs.by.peaks <- with(converted, split(segments, segments$peaks))
  for(peaks in names(segs.by.peaks)){
    model.segs <- segs.by.peaks[[peaks]]
    if(is.data.frame(model.segs)){
      segs.by.peaks.fit[[peaks]][[fit.name]] <-
        data.frame(fit.name, model.segs)
    }
  }
}
do.call(rbind, loss.list)

segs1 <- do.call(rbind, segs.by.peaks.fit[["1"]])
breaks1 <- subset(segs1, min(chromStart) < chromStart)
library(ggplot2)
ggplot()+
  ggtitle(paste("PeakSegJointSeveral runs PeakSegJointHeuristic",
                "and keeps only the most likely model"))+
  geom_step(aes(chromStart/1e3, count),
            color="grey50",
            data=H3K4me3.TDH.other.chunk8)+
  geom_vline(aes(xintercept=chromStart/1e3),
             data=breaks1,
             color="green",
             linetype="dashed")+
  geom_segment(aes(chromStart/1e3, mean,
                   xend=chromEnd/1e3, yend=mean),
               size=1,
               color="green",
               data=segs1)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(sample.id ~ fit.name, labeller=function(var, val){
    if(var=="sample.id"){
      sub("McGill0", "", val)
    }else{
      ifelse(val=="several", "best of several", paste("bin.factor =", val))
    }
  }, scales="free")

segs.by.peaks <- list()
for(peaks in 8:10){
  segs.by.peaks[[paste(peaks)]] <-
    data.frame(peaks, segs.by.peaks.fit[[paste(peaks)]][["several"]])
}
segs <- do.call(rbind, segs.by.peaks)
breaks <- subset(segs, min(chromStart) < chromStart)
library(ggplot2)
ggplot()+
  ggtitle("PeakSegJoint models with 8-10 peaks")+
  geom_step(aes(chromStart/1e3, count),
            color="grey50",
            data=H3K4me3.TDH.other.chunk8)+
  geom_vline(aes(xintercept=chromStart/1e3),
             data=breaks,
             color="green",
             linetype="dashed")+
  geom_segment(aes(chromStart/1e3, mean,
                   xend=chromEnd/1e3, yend=mean),
               size=1,
               color="green",
               data=segs)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(sample.id ~ peaks, labeller=function(var, val){
    if(var=="sample.id"){
      sub("McGill0", "", val)
    }else{
      paste(val, "peaks")
    }
  }, scales="free")

}
