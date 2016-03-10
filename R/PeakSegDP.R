PeakSegDP <- structure(function
### Compute the PeakSeg model on a data.frame of compressed sequence
### reads.
(compressed,
### data.frame with columns chromStart, chromEnd, count.
 maxPeaks
### maximum number of peaks to consider.
 ){
  stopifnot(diff(compressed$chromStart) > 0)
  count <- compressed$count
  weight <- compressed$bases <- with(compressed, chromEnd-chromStart)
  stopifnot(is.integer(count))
  stopifnot(is.integer(weight))
  stopifnot(count >= 0)
  stopifnot(weight > 0)
  stopifnot(is.integer(maxPeaks))
  stopifnot(length(maxPeaks) == 1)
  stopifnot(length(count) == length(weight))
  maxSegments <- maxPeaks * 2 + 1
  stopifnot(maxSegments > 0)
  stopifnot(maxSegments <= nrow(compressed))
  fit <- cDPA(count, weight, maxSegments)
  segment.ends <- getPath(fit)
  results <- list(breaks=list())
  for(peaks in 0:maxPeaks){
    peak.list <- list()
    segments <- as.integer(peaks*2 + 1)
    model.i <- peaks * 2 + 1
    last.i <- as.integer(segment.ends[model.i, 1:model.i])
    not.feasible <- any(last.i < 1 | is.na(last.i))
    if(not.feasible){
      s <- ifelse(peaks==1, "", "s")
      warning("infeasible model with ", peaks, " peak", s)
    }else{
      break.after <- last.i[-model.i]
      first.i <- as.integer(c(1, break.after+1))
      model.error <- 0
      if(length(break.after)){
        results$breaks[[paste(peaks)]] <-
          data.frame(peaks, segments, break.after,
                     chromEnd=compressed$chromEnd[break.after])
      }
      for(segment.i in seq_along(last.i)){
        status <- ifelse(segment.i %% 2, "background", "peak")
        first <- first.i[[segment.i]]
        last <- last.i[[segment.i]]
        seg.data <- compressed[first:last,]
        count.num <- as.numeric(seg.data$count)
        bases.num <- as.numeric(seg.data$bases)
        seg.mean <- sum(count.num*bases.num)/sum(bases.num)
        model.error <- model.error + with(seg.data, {
          PoissonLoss(count, seg.mean, bases)
        })
        chromStart <- seg.data$chromStart[1]
        chromEnd <- seg.data$chromEnd[nrow(seg.data)]
        results$segments[[paste(peaks, segment.i)]] <- 
          data.frame(mean=seg.mean,
                     first,
                     last,
                     chromStart,
                     chromEnd,
                     status,
                     peaks,
                     segments)
        if(status == "peak"){
          peak.list[[paste(segment.i)]] <-
            data.frame(first, last,
                       chromStart, chromEnd,
                       peaks, segments)
        }
      }#segment.i
      results$peaks[[as.character(peaks)]] <- do.call(rbind, peak.list)
      results$error[[as.character(peaks)]] <- 
        data.frame(segments=model.i, peaks, error=model.error)
    }#!not.feasible
  }#peaks
  no.peaks <-
    data.frame(first=integer(), last=integer(),
               chromStart=integer(), chromEnd=integer(),
               peaks=integer(), segments=integer())
  results$peaks <- c(list("0"=no.peaks), results$peaks)
  results$error <- do.call(rbind, results$error)
  results$segments <- do.call(rbind, results$segments)
  results$breaks <- do.call(rbind, results$breaks)
  results
}, ex=function(){
  data(chr11ChIPseq)
  one <- subset(chr11ChIPseq$coverage, sample.id=="McGill0002")
  fit <- PeakSegDP(one, 3L)
  library(ggplot2)
  ggplot()+
    geom_step(aes(chromStart/1e3, count), data=one)+
    geom_segment(aes(chromStart/1e3, mean,
                     xend=chromEnd/1e3, yend=mean),
                 data=fit$segments, color="green")+
    geom_segment(aes(chromStart/1e3, 0,
                     xend=chromEnd/1e3, yend=0),
                 data=subset(fit$segments, status=="peak"),
                 size=3, color="deepskyblue")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(peaks ~ ., scales="free", labeller=function(var, val){
      s <- ifelse(val==1, "", "s")
      paste0(val, " peak", s)
    })
})
