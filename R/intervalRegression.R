exactModelSelection <- structure(function # Exact model selection function
### Given a set of optimal costs C_i, and model complexity values K_i,
### and a model selection function i*(lambda) = argmin_i C_i +
### lambda*K_i, compute a set of consecutive (K_i, min.lambda,
### max.lambda) with i being the solution for every lambda in
### (min.lambda, max.lambda).
(cost,
### numeric vector: optimal costs C_i.
 model.complexity,
### numeric vector: model complexity K_i.
 peaks){
  stopifnot(is.numeric(cost))
  stopifnot(is.numeric(model.complexity))
  stopifnot(diff(model.complexity) > 0)
  stopifnot(diff(cost) < 0)
  stopifnot(length(cost) == length(model.complexity))
  n.models <- length(cost)
  Kmax <- model.complexity[n.models]
  Kcurrent <- Kmax
  Lcurrent <- 0
  vK <- Kmax
  vL <- 0
  vP <- peaks[n.models]
  i <- 2
  min.complexity <- model.complexity[1]
  while(Kcurrent > min.complexity) {
    is.smaller <- model.complexity < Kcurrent
    is.current <- model.complexity == Kcurrent
    smallerK <- model.complexity[is.smaller]
    smallerPeaks <- peaks[is.smaller]
    cost.term <- cost[is.current] - cost[is.smaller]
    complexity.term <- smallerK - model.complexity[is.current]
    lambdaTransition <- cost.term/complexity.term
    next.i <- which.min(lambdaTransition)
    Kcurrent <- smallerK[next.i]
    Lcurrent <- min(lambdaTransition)
    vL[i] <- Lcurrent
    vK[i] <- Kcurrent
    vP[i] <- smallerPeaks[next.i]
    i <- i + 1
  }
  L <- log(vL)
  data.frame(min.log.lambda = L,
             max.log.lambda = c(L[-1], Inf),
             model.complexity = vK,
             peaks=vP,
             min.lambda = vL,
             max.lambda = c(vL[-1], Inf))
},ex=function(){
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
                  label=sprintf("%d peak%s optimal", peaks,
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
})

largestContinuousMinimum <- structure(function
### Find the run of minimum cost with the largest size.
(cost,
 size
 ){
  m <- min(cost)
  is.min <- cost == m
  d <- c(diff(c(FALSE,is.min,FALSE)))
  ##print(data.frame(cost=c(cost,NA),size=c(size,NA),diff=d))
  starts <- which(d==1)
  ends <- which(d==-1)-1
  runs <- data.frame(starts,ends)
  ##print(runs)
  runs$size <- sapply(seq_along(starts),function(i){
    sum(size[ starts[i]:ends[i] ])
  })
  ##print(runs)
  if(1 < sum(runs$size==Inf)){
    list(start=1, end=length(cost))
  }else{
    largest <- which.max(runs$size)
    list(start=starts[largest],end=ends[largest])
  }
}, ex=function(){
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
    several.samples <- if(peaks.str %in% names(peaks.by.peaks)){
      peaks.by.peaks[[peaks.str]]
    }else{
      Peaks()
    }
    peaks.by.sample <- split(several.samples, several.samples$sample.id)
    error.by.sample <- list()
    for(sample.id in names(regions.by.sample)){
      one.sample.peaks <- if(sample.id %in% names(peaks.by.sample)){
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
})

