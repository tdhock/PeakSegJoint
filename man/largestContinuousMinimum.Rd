\name{largestContinuousMinimum}
\alias{largestContinuousMinimum}
\title{largestContinuousMinimum}
\description{Find the run of minimum cost with the largest size.}
\usage{largestContinuousMinimum(cost, size)}
\arguments{
  \item{cost}{
}
  \item{size}{
}
}



\author{Toby Dylan Hocking}




\examples{
data(H3K36me3.TDH.other.chunk1)
lims <- c(43000000, 43200000) # left
some.counts <-
  subset(H3K36me3.TDH.other.chunk1$counts,
         lims[1] < chromEnd & chromStart < lims[2])
fit <- PeakSegJointHeuristic(some.counts)
converted <- ConvertModelList(fit)
## Ignore a model if there is another one with lower peaks and loss.
all.loss <- converted$loss
all.loss$cummin <- cummin(all.loss$loss)
some.loss <- subset(all.loss, loss == cummin)
## Compute PeakError on this sequence of models.
some.regions <- 
  subset(H3K36me3.TDH.other.chunk1$regions,
         chromStart < lims[2])
regions.by.sample <- split(some.regions, some.regions$sample.id)
peaks.by.peaks <- split(converted$peaks, converted$peaks$peaks)
library(PeakError)
error.list <- list()
for(peaks in some.loss$peaks){
  peaks.str <- paste(peaks)
  several.samples <- if(peaks.str \%in\% names(peaks.by.peaks)){
    peaks.by.peaks[[peaks.str]]
  }else{
    Peaks()
  }
  peaks.by.sample <- split(several.samples, several.samples$sample.id)
  error.by.sample <- list()
  for(sample.id in names(regions.by.sample)){
    one.sample.peaks <- if(sample.id \%in\% names(peaks.by.sample)){
      peaks.by.sample[[sample.id]]
    }else{
      Peaks()
    }
    one.sample.regions <- regions.by.sample[[sample.id]]
    error <- PeakErrorChrom(one.sample.peaks, one.sample.regions)
    error.by.sample[[sample.id]] <- data.frame(sample.id, error)
  }
  peaks.error <- do.call(rbind, error.by.sample)
  fp <- sum(peaks.error$fp)
  fn <- sum(peaks.error$fn)
  error.list[[peaks.str]] <-
    data.frame(peaks, errors=fp+fn, regions=nrow(peaks.error))
}
error.sum <- do.call(rbind, error.list)
show.error <- data.frame(what="incorrect regions", error.sum)
show.loss <- data.frame(what="Poisson loss", some.loss)
library(ggplot2)
ggplot()+
  geom_point(aes(peaks, loss),
             data=show.loss)+
  geom_line(aes(peaks, loss),
            data=show.loss)+
  geom_point(aes(peaks, errors),
             data=show.error)+
  geom_line(aes(peaks, errors),
            data=show.error)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(what ~ ., scales="free")
## Compute the exact model selection function.
exact.df <- with(some.loss, exactModelSelection(loss, peaks, peaks))
## Compute the target interval.
rownames(show.error) <- show.error$peaks
exact.df$errors <- show.error[paste(exact.df$peaks), "errors"]
indices <- with(exact.df, {
  largestContinuousMinimum(errors, max.log.lambda-min.log.lambda)
})
target.interval <-
  data.frame(min.log.lambda=exact.df$min.log.lambda[indices$start],
             max.log.lambda=exact.df$max.log.lambda[indices$end])
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(what ~ ., scales="free")+
  geom_tallrect(aes(xmin=min.log.lambda, xmax=max.log.lambda),
                data=target.interval,
                fill="grey",
                alpha=0.5)+
  ggtitle(paste("target interval of penalty values with",
                "minimal incorrect regions in grey"))+
  geom_segment(aes(min.log.lambda, errors,
                   xend=max.log.lambda, yend=errors),
               data=data.frame(exact.df, what="incorrect regions"))+
  geom_segment(aes(min.log.lambda, peaks,
                   xend=max.log.lambda, yend=peaks),
               data=data.frame(exact.df, what="peaks"))
}
