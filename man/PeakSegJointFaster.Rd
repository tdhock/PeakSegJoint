\name{PeakSegJointFaster}
\alias{PeakSegJointFaster}
\title{PeakSegJointFaster}
\description{Run the PeakSegJointFaster heuristic optimization algorithm, for
several bin.factor parameter values, keeping only the most likely
model found. This gives an approximate solution to a multi-sample
Poisson maximum likelihood segmentation problem. Given S samples,
this function computes a sequence of S+1 PeakSegJoint models, with
0, ..., S samples with an overlapping peak (maximum of one peak
per sample). It also computes for G groups, the seq of G+1
models, with 0, ..., G groups with an overlapping peak.}
\usage{PeakSegJointFaster(profiles, 
    bin.factor.vec = 2:7)}
\arguments{
  \item{profiles}{data.frame with columns sample.id, sample.group, chromStart,
chromEnd, count.}
  \item{bin.factor.vec}{Size of bin pyramid. Bigger values result in slower computation.}
}

\value{List of model fit results.}

\author{Toby Dylan Hocking}




\examples{

library(PeakSegJoint)
data(H3K36me3.TDH.other.chunk1, envir=environment())
some.counts <- subset(
  H3K36me3.TDH.other.chunk1$counts,
  43000000 < chromEnd &
  chromStart < 43200000)
some.counts$sample.group <- some.counts$cell.type

fit <- PeakSegJointFaster(some.counts, 2:7)

if(interactive() && require(ggplot2)){

  both <- with(fit, rbind(
    data.frame(model="sample", sample.modelSelection),
  data.frame(model="group", group.modelSelection)))
  ggplot()+
    ggtitle("model selection functions")+
    scale_size_manual(values=c(sample=2, group=1))+
    geom_segment(aes(min.log.lambda, complexity,
                     color=model, size=model,
                     xend=max.log.lambda, yend=complexity),
                 data=both)+
    xlab("log(penalty)")+
    ylab("model complexity (samples or groups with a common peak)")

}

}
