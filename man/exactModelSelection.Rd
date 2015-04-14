\name{exactModelSelection}
\alias{exactModelSelection}
\title{Exact model selection function}
\description{Given a set of optimal costs C_i, and model complexity values K_i,
and a model selection function i*(lambda) = argmin_i C_i +
lambda*K_i, compute a set of consecutive (K_i, min.lambda,
max.lambda) with i being the solution for every lambda in
(min.lambda, max.lambda).}
\usage{exactModelSelection(cost, model.complexity, peaks)}
\arguments{
  \item{cost}{numeric vector: optimal costs C_i.}
  \item{model.complexity}{numeric vector: model complexity K_i.}
  \item{peaks}{
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
## Calculate the exact path of breakpoints in the optimal number of
## peaks function.
exact.df <- with(some.loss, exactModelSelection(loss, peaks, peaks))
rownames(some.loss) <- some.loss$peaks
intercept <- some.loss[as.character(exact.df$peaks), "loss"]
exact.df$cost <- intercept + exact.df$min.lambda * exact.df$model.complexity
exact.df$next.cost <- c(exact.df$cost[-1], NA)
library(ggplot2)
ggplot()+
  geom_point(aes(min.lambda, cost),
             data=exact.df, pch=1, color="red")+
  geom_segment(aes(min.lambda, cost,
                   xend=max.lambda, yend=next.cost),
               data=exact.df, color="red", size=1.5)+
  geom_text(aes((min.lambda+max.lambda)/2, (cost+next.cost)/2,
                label=sprintf("\%d peak\%s optimal", peaks,
                  ifelse(peaks==1, "", "s"))),
            data=exact.df, color="red", hjust=0, vjust=1.5)+
  geom_abline(aes(slope=peaks, intercept=loss), data=some.loss)+
  geom_text(aes(0, loss, label=peaks),
            data=some.loss, hjust=1.5, color="red")+
  ggtitle("model selection: cost = loss_k + lambda*segments_k")
## Solve the optimization using grid search.
L.grid <- with(exact.df,{
  seq(min(max.log.lambda)-1,
      max(min.log.lambda)+1,
      l=100)
})
lambda.grid <- exp(L.grid)
kstar.grid <- sapply(lambda.grid,function(lambda){
  crit <- some.loss$peaks * lambda + some.loss$loss
  picked <- which.min(crit)
  some.loss$peaks[picked]
})
grid.df <- data.frame(log.lambda=L.grid, peaks=kstar.grid)
## Compare the results.
ggplot()+
  geom_segment(aes(min.log.lambda, peaks,
                   xend=max.log.lambda, yend=peaks),
               data=exact.df)+
  geom_point(aes(log.lambda, peaks),
             data=grid.df, color="red", pch=1)+
  ylab("optimal model complexity (peaks)")+
  xlab("log(lambda)")
}
