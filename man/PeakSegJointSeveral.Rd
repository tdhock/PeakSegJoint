\name{PeakSegJointSeveral}
\alias{PeakSegJointSeveral}
\title{PeakSegJointSeveral}
\description{Run PeakSegJoint with several different bin.factor values, keeping
only the models with lowest Poisson loss for each peak size.}
\usage{PeakSegJointSeveral(profiles, bin.factors = c(2, 3, 5))}
\arguments{
  \item{profiles}{data.frame or list of them from ProfileList.}
  \item{bin.factors}{integer vector of suboptimality parameters, bigger values are
slower.}
}

\value{Model fit list as is returned by PeakSegJointHeuristic.}

\author{Toby Dylan Hocking}




\examples{
data(H3K4me3.TDH.other.chunk8)
bf.vec <- c(2, 3, 5)
fit.list <-
  list(several=PeakSegJointSeveral(H3K4me3.TDH.other.chunk8, bf.vec))
for(bf in bf.vec){
  fit.list[[paste(bf)]] <-
    PeakSegJointHeuristic(H3K4me3.TDH.other.chunk8, bf)
}
loss.list <- list()
segs.list <- list()
for(fit.name in names(fit.list)){
  fit <- fit.list[[fit.name]]
  loss.list[[fit.name]] <- sapply(fit$models, "[[", "loss")
  converted <- ConvertModelList(fit)
  segs.list[[fit.name]] <-
    data.frame(fit.name,
               subset(converted$segments, peaks == 2))
}
do.call(rbind, loss.list)
segs <- do.call(rbind, segs.list)
breaks <- subset(segs, min(chromStart) < chromStart)
library(ggplot2)
ggplot()+
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
  facet_grid(sample.id ~ fit.name, labeller=function(var, val){
    sub("McGill0", "", val)
  }, scales="free")
}
