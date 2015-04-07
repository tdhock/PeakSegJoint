multiSampleSegHeuristic <- structure(function
### Find one peak common to several samples.
(profiles,
### List of data.frames with columns chromStart, chromEnd, count, or
### single data.frame with additional column sample.id.
 bin.factor=2L
### Size of bin pyramid.
 ){
  stopifnot(is.numeric(bin.factor))
  stopifnot(length(bin.factor)==1)
  if(is.data.frame(profiles)){
    stopifnot(!is.null(profiles$sample.id))
    profiles <- split(profiles, profiles$sample.id, drop=TRUE)
  }
  stopifnot(is.list(profiles))
  for(profile.i in seq_along(profiles)){
    df <- profiles[[profile.i]]
    stopifnot(is.data.frame(df))
    stopifnot(is.integer(df$chromStart))
    stopifnot(is.integer(df$chromEnd))
    stopifnot(is.integer(df$count))
    profiles[[profile.i]] <- df[, c("chromStart", "chromEnd", "count")]
  }
  chromStartEnd <-
    .Call("multiSampleSegHeuristic_interface",
          profiles,
          as.integer(bin.factor),
          PACKAGE="PeakSegDP")
  data.frame(chromStart=chromStartEnd[1],
             chromEnd=chromStartEnd[2])
}, ex=function(){
  library(PeakSegDP)

  data(chr11ChIPseq)
  two <- subset(chr11ChIPseq$coverage,
                118090000 < chromStart &
                chromEnd < 118100000 &
                sample.id %in% c("McGill0002", "McGill0004"))
  ## Find the best peak location across 2 samples.
  optimal.seconds <- system.time({
    optimal <- multiSampleSegOptimal(two)
  })[["elapsed"]]
  heuristic.seconds <- system.time({
    heuristic <- multiSampleSegHeuristic(two, 2)
  })[["elapsed"]]
  rbind(heuristic.seconds, optimal.seconds)
  peaks <-
    rbind(data.frame(optimal, model="optimal"),
          data.frame(heuristic, model="heuristic"))
  library(ggplot2)
  ggplot()+
    geom_step(aes(chromStart/1e3, count), data=two)+
    scale_size_manual(values=c(optimal=2, heuristic=1))+
    geom_segment(aes(chromStart/1e3, 0,
                     color=model, size=model,
                     xend=chromEnd/1e3, yend=0),
                 data=peaks)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., scales="free")

  ## Find the best peak location across 4 samples.
  four <- subset(chr11ChIPseq$coverage,
                 118120000 < chromStart &
                 chromEnd < 118126000) 
  optimal.seconds <- system.time({
    optimal <- multiSampleSegOptimal(four)
  })[["elapsed"]]
  heuristic.seconds <- system.time({
    heuristic <- multiSampleSegHeuristic(four, 2)
  })[["elapsed"]]
  rbind(heuristic.seconds, optimal.seconds)
  peaks <-
    rbind(data.frame(optimal, model="optimal"),
          data.frame(heuristic, model="heuristic"))
  ggplot()+
    geom_step(aes(chromStart/1e3, count), data=four)+
    scale_size_manual(values=c(optimal=2, heuristic=1))+
    geom_segment(aes(chromStart/1e3, 0,
                     color=model, size=model,
                     xend=chromEnd/1e3, yend=0),
                 data=peaks)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., scales="free")

  ## A fake data set with two profiles with very different scales for
  ## the count variable, showing that the profile with the small
  ## counts will be basically ignored when computing the optimal
  ## peak.
  count <-
    c(rep(c(0, 1), each=100),
      rep(c(10, 11), each=200),
      rep(c(0, 1), each=100))
  chromEnd <- seq_along(count)
  chromStart <- chromEnd-1L
  offset <- 50L
  multi.scale <- 
  rbind(data.frame(sample.id="low", chromStart, chromEnd,
                   count=as.integer(count)),
        data.frame(sample.id="hi",
                   chromStart=chromStart+offset,
                   chromEnd=chromEnd+offset,
                   count=as.integer(count*1000)))
  heuristic.seconds <- system.time({
    heuristic <- multiSampleSegHeuristic(multi.scale, 2)
  })[["elapsed"]]
  optimal.seconds <- system.time({
    optimal <- multiSampleSegOptimal(multi.scale)
  })[["elapsed"]]
  rbind(heuristic.seconds, optimal.seconds)
  peaks <-
    rbind(data.frame(optimal, model="optimal"),
          data.frame(heuristic, model="heuristic"))
  ggplot()+
    geom_step(aes(chromStart, count), data=multi.scale)+
    scale_size_manual(values=c(optimal=2, heuristic=1))+
    geom_segment(aes(chromStart, 0,
                     color=model, size=model,
                     xend=chromEnd, yend=0),
                 data=peaks)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., scales="free")
  
  ## plot microbenchmark time versus parameter.
  data(H3K4me3.TDH.immune.chunk12.cluster4)
  many <- H3K4me3.TDH.immune.chunk12.cluster4
  optimal <- c(27998215, 27999159)
  optimal.seconds <- system.time({
    optimal <- multiSampleSegOptimal(many)
  })[["elapsed"]]
  library(microbenchmark)
  m.args <- list()
  results <- list()
  param.values <-
    c(2, 3, 5, 10, 15, 20, 25, 30, 50, 75, 100, 125, 150,
      200, 300, 400, 500, 600)
  for(param in param.values){
    param.name <- paste(param)
    param.int <- as.integer(param)
    m.args[[param.name]] <- substitute({
      print(PINT)
      peak <- multiSampleSegHeuristic(many, PINT)
      diff.bases <- sum(abs(peak - optimal))
      results[[paste(PINT)]] <- data.frame(param=PINT, peak, diff.bases)
    }, list(PINT=param.int))
  }
  times <- microbenchmark(list=m.args, times=3)
  results.df <- do.call(rbind, results)
  optimal.timing <- 
    data.frame(seconds=optimal.seconds, what="seconds")
  ggplot()+
    ggtitle(paste("heuristic segmentation accurate within 10 bases",
                  "and much faster than optimal segmentation",
                  sep="\n"))+
    scale_x_log10()+
    scale_y_log10()+
    geom_point(aes(as.numeric(as.character(param)), diff.bases),
               data=data.frame(results.df, what="diff.bases"),
               pch=1)+
    facet_grid(what ~ ., scales="free")+
    geom_hline(aes(yintercept=seconds), data=optimal.timing)+
    geom_text(aes(10, seconds,
                  label=sprintf("optimal segmentation = %.1f seconds",
                    seconds)),
              hjust=0,
              vjust=1.5,
              data=optimal.timing)+
    geom_point(aes(as.numeric(as.character(expr)), time/1e9),
               data=data.frame(times, what="seconds"), pch=1)
  ## The biggest observed deviation of the Heuristic result from the
  ## optimal result is 17 bases (for bin.factor=30).
  stopifnot(results.df$diff.bases < 20)

  heuristic <- subset(results.df, param=="2", select=c(chromStart, chromEnd))
  peaks <-
    rbind(data.frame(optimal, model="optimal"),
          data.frame(heuristic, model="heuristic"))
  ggplot()+
    geom_step(aes(chromStart/1e3, count), data=many)+
    scale_size_manual(values=c(optimal=2, heuristic=1))+
    geom_segment(aes(chromStart/1e3, 0,
                     color=model, size=model,
                     xend=chromEnd/1e3, yend=0),
                 data=peaks)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., scales="free")

  ## Heuristic in R code for debugging.
  profiles <- many
  profile.list <- split(profiles, profiles$sample.id, drop=TRUE)
  max.chromStart <- max(sapply(profile.list, with, chromStart[1]))
  min.chromEnd <- min(sapply(profile.list, with, chromEnd[length(chromEnd)]))
  c(max.chromStart, min.chromEnd)
  bases <- min.chromEnd-max.chromStart
  n.samples <- length(profile.list)
  small.chromEnd <- (max.chromStart+1):min.chromEnd
  small.bins <-
    matrix(NA, bases, n.samples,
           dimnames=list(position=small.chromEnd,
             sample.id=names(profile.list)))
  for(sample.id in names(profile.list)){
    one <- profile.list[[sample.id]]
    bins <- binSum(one, max.chromStart, n.bins=bases)
    stopifnot(bins$chromEnd == small.chromEnd)
    small.bins[, sample.id] <- bins$count
  }
  bin.factor <- 10L
  bases.per.bin <- 1L
  while(bases/bases.per.bin/bin.factor >= 4){
    bases.per.bin <- bases.per.bin * bin.factor
  }
  n.bins <- as.integer(bases %/% bases.per.bin + 1L)
  na.mat <- 
    matrix(NA, n.bins, n.samples,
           dimnames=list(bin=NULL, sample.id=names(profile.list)))
  first.cumsums <- list(count=na.mat)
  bin.list <- list()
  norm.list <- list()
  for(sample.i in seq_along(profile.list)){
    sample.id <- names(profile.list)[sample.i]
    one <- profile.list[[sample.i]]
    max.count <- max(one$count)
    bins <- binSum(one, max.chromStart, bases.per.bin, n.bins)
    stopifnot(n.bins == nrow(bins))
    bins[nrow(bins), "chromEnd"] <- min.chromEnd
    bins$mean <- with(bins, count/(chromEnd-chromStart))
    bins$mean.norm <- bins$mean/max.count
    bin.list[[sample.id]] <- data.frame(sample.id, rbind(bins, NA))
    bases.vec <- with(bins, chromEnd-chromStart)
    stopifnot(bins$count >= 0)
    first.cumsums$count[, sample.i] <- cumsum(bins$count)
    one$count.norm <- one$count/max.count
    norm.list[[sample.i]] <- one
  }
  bin.df <- do.call(rbind, bin.list)
  norm.df <- do.call(rbind, norm.list)
  ## The formula for the optimal Poisson loss 
  ## for 1 segment with d integer data points x_j is
  ## \sum_{j=1}^d m - x_j \log m_j =
  ##   ( \sum_{j=1}^d x_j ) (1-\log m)
  ## where the segment mean m = (\sum x_j)/d,
  OptimalPoissonLoss <- function(mean.value, cumsum.value){
    ifelse(mean.value == 0, 0, cumsum.value * (1-log(mean.value)))
  }
  loss.list <- list()
  seg.list <- list()
  for(seg1.last in 1:(n.bins-2)){
    seg1.cumsums <- first.cumsums$count[seg1.last, ]
    seg1.bases <- seg1.last*bases.per.bin
    seg1.chromEnd <- seg1.bases + max.chromStart
    seg1.means <- seg1.cumsums/seg1.bases
    seg1.loss.vec <- OptimalPoissonLoss(seg1.means, seg1.cumsums)
    seg1.loss <- sum(seg1.loss.vec)
    for(seg2.last in (seg1.last+1):(n.bins-1)){
      cumsum.seg2.end <- first.cumsums$count[seg2.last, ]
      seg2.cumsums <- cumsum.seg2.end-seg1.cumsums
      seg12.bases <- seg2.last*bases.per.bin
      seg2.bases <- seg12.bases - seg1.bases
      seg2.chromEnd <- seg1.chromEnd + seg2.bases
      seg2.means <- seg2.cumsums/seg2.bases
      seg2.loss.vec <- OptimalPoissonLoss(seg2.means, seg2.cumsums)
      seg2.loss <- sum(seg2.loss.vec)
      
      seg3.cumsums <- first.cumsums$count[n.bins, ]-cumsum.seg2.end
      seg3.bases <- bases - seg12.bases
      seg3.means <- seg3.cumsums/seg3.bases
      seg3.loss.vec <- OptimalPoissonLoss(seg3.means, seg3.cumsums)
      seg3.loss <- sum(seg3.loss.vec)
      total.loss <- seg1.loss + seg2.loss + seg3.loss
      loss.list[[paste(seg1.last, seg2.last)]] <-
        data.frame(seg1.last, seg2.last, total.loss,
                   seg1.loss, seg2.loss, seg3.loss,
                   peakStart=seg1.last*bases.per.bin+max.chromStart,
                   peakEnd=seg2.last*bases.per.bin+max.chromStart)
      mean.mat <- rbind(seg1.means, seg2.means, seg3.means)
      for(sample.id in colnames(mean.mat)){
        seg.list[[paste(seg1.last, seg2.last)]][[sample.id]] <-
          data.frame(sample.id,
                     chromStart=c(max.chromStart, seg1.chromEnd, seg2.chromEnd),
                     chromEnd=c(seg1.chromEnd, seg2.chromEnd, min.chromEnd),
                     mean=mean.mat[, sample.id])
      }
    }
  }
  loss.df <- do.call(rbind, loss.list)
  loss.best <- loss.df[which.min(loss.df$total.loss), ]
  seg.df <- do.call(rbind, {
    seg.list[[with(loss.best, paste(seg1.last, seg2.last))]]
  })
  
  ggplot()+
    scale_color_manual(values=c(data="grey50", bins="black", segments="green"))+
    geom_step(aes(chromStart/1e3, count, color=what),
              data=data.frame(norm.df, what="data"))+
    geom_segment(aes(chromStart/1e3, mean,
                     xend=chromEnd/1e3, yend=mean,
                     color=what),
                 size=1,
                 data=data.frame(seg.df, what="segments"))+
    geom_segment(aes(chromStart/1e3, mean,
                     xend=chromEnd/1e3, yend=mean,
                     color=what),
                 data=data.frame(bin.df, what="bins"))+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
      sub("McGill0", "", val)
    })

  ## These need to be computed once and then they never change.
  last.cumsums <- list()
  before.cumsums <- list()
  for(data.type in names(first.cumsums)){
    mat.or.vec <- first.cumsums[[data.type]]
    if(is.matrix(mat.or.vec)){
      last.cumsums[[data.type]] <- mat.or.vec[nrow(mat.or.vec), ]
      before.cumsums$left[[data.type]] <- if(loss.best$seg1.last == 1){
        rep(0, n.samples)
      }else{
        mat.or.vec[loss.best$seg1.last-1, ]
      }
      before.cumsums$right[[data.type]] <- mat.or.vec[loss.best$seg2.last-1, ]
    }else{
      last.cumsums[[data.type]] <- mat.or.vec[length(mat.or.vec)]
      before.cumsums$left[[data.type]] <- if(loss.best$seg1.last == 1){
        0
      }else{
        mat.or.vec[loss.best$seg1.last-1]
      }
      before.cumsums$right[[data.type]] <- mat.or.vec[loss.best$seg2.last-1]
    }
  }
  last.chromEnd <- min.chromEnd
  n.bins.zoom <- bin.factor * 2L
  n.cumsum <- n.bins.zoom + 1L

  ## These will change at the end of each iteration.
  peakStart <- loss.best$peakStart
  peakEnd <- loss.best$peakEnd
  search.list <- list()
  best.list <- list()
  while(bases.per.bin > 1){
    left.chromStart <- peakStart - bases.per.bin
    right.chromStart <- peakEnd-bases.per.bin
    bases.per.bin <- as.integer(bases.per.bin / bin.factor)

    cumsum.list <- function(chromStart){
      limits <- bases.per.bin*(0:(n.cumsum-1))+chromStart
      chromStart.vec <- limits[-length(limits)]
      chromEnd.vec <- limits[-1]
      intervals <-
        paste0(chromStart.vec, "-", chromEnd.vec)
      dn <- list(bin=c("before", intervals), sample.id=names(profile.list))
      m <- matrix(NA, n.cumsum, n.samples, dimnames=dn)
      list(count=m, 
           chromStart=chromStart.vec, chromEnd=chromEnd.vec)
    }
    cumsum.mats <-
      list(left=cumsum.list(left.chromStart),
           right=cumsum.list(right.chromStart))

    for(sample.i in seq_along(profile.list)){
      one <- profile.list[[sample.i]]
      lr.list <-
        list(left=binSum(one, left.chromStart, bases.per.bin, n.bins.zoom),
             right=binSum(one, right.chromStart, bases.per.bin, n.bins.zoom,
                          empty.as.zero=TRUE))
      for(lr in names(lr.list)){
        lr.bins <- lr.list[[lr]]
        stopifnot(nrow(lr.bins) == n.bins.zoom)
        lr.bases <- with(lr.bins, chromEnd-chromStart)
        lr.before <- before.cumsums[[lr]]
        lr.counts <- list(count=lr.bins$count)
        for(data.type in names(lr.counts)){
          lr.count.vec <-
            c(lr.before[[data.type]][[sample.i]], lr.counts[[data.type]])
          cumsum.mats[[lr]][[data.type]][, sample.i] <-
            cumsum(lr.count.vec)
        }
      }
    }
    possible.grid <- 
      expand.grid(left.cumsum.row=3:n.cumsum, right.cumsum.row=2:n.cumsum)
    possible.grid$left.chromStart <-
      cumsum.mats$left$chromStart[possible.grid$left.cumsum.row-1]
    possible.grid$left.chromEnd <-
      cumsum.mats$left$chromEnd[possible.grid$left.cumsum.row-1]
    possible.grid$right.chromStart <-
      cumsum.mats$right$chromStart[possible.grid$right.cumsum.row-1]
    possible.grid$right.chromEnd <-
      cumsum.mats$right$chromEnd[possible.grid$right.cumsum.row-1]
    feasible.grid <-
      subset(possible.grid,
             left.chromEnd <= right.chromStart &
             right.chromEnd < last.chromEnd)
    feasible.grid$model.i <- 1:nrow(feasible.grid)
    model.list <- list()
    seg.list <- list()
    for(model.i in feasible.grid$model.i){
      model.row <- feasible.grid[model.i, ]

      seg1.i <- model.row$left.cumsum.row-1
      seg1.cumsums <- cumsum.mats$left$count[seg1.i, , drop=FALSE]
      seg1.chromEnd <- cumsum.mats$left$chromStart[seg1.i]
      seg1.bases <- seg1.chromEnd-max.chromStart
      seg1.means <- seg1.cumsums/seg1.bases
      seg1.loss <- OptimalPoissonLoss(seg1.means, seg1.cumsums)
      seg.list[[paste(model.i, 1)]] <-
        data.frame(chromStart=max.chromStart, chromEnd=seg1.chromEnd,
                   mean=seg1.means, sample.id=names(profile.list),
                   model.i,
                   row.names=NULL)

      seg1.mat <-
        small.bins[small.chromEnd <= seg1.chromEnd, , drop=FALSE]
      stopifnot(nrow(seg1.mat) == seg1.bases)
      stopifnot(all.equal(as.numeric(colMeans(seg1.mat)),
                          as.numeric(seg1.means)))

      cumsum.seg2.end <-
        cumsum.mats$right$count[model.row$right.cumsum.row, , drop=FALSE]
      seg2.cumsums <- cumsum.seg2.end-seg1.cumsums
      seg2.chromEnd <- cumsum.mats$right$chromEnd[model.row$right.cumsum.row-1]
      seg2.bases <- seg2.chromEnd-seg1.chromEnd
      seg2.means <- seg2.cumsums/seg2.bases
      seg2.loss <- OptimalPoissonLoss(seg2.means, seg2.cumsums)
      seg.list[[paste(model.i, 2)]] <-
        data.frame(chromStart=seg1.chromEnd, chromEnd=seg2.chromEnd,
                   mean=seg2.means, sample.id=names(profile.list),
                   model.i,
                   row.names=NULL)

      is.seg2 <-
        seg1.chromEnd < small.chromEnd &
          small.chromEnd <= seg2.chromEnd
      stopifnot(sum(is.seg2) == seg2.bases)
      seg2.mat <- small.bins[is.seg2, , drop=FALSE]
      stopifnot(all.equal(as.numeric(colMeans(seg2.mat)),
                          as.numeric(seg2.means)))
      
      seg3.cumsums <- last.cumsums$count-cumsum.seg2.end
      seg3.bases <- last.chromEnd-seg2.chromEnd
      seg3.means <- seg3.cumsums/seg3.bases
      seg3.loss <- OptimalPoissonLoss(seg3.means, seg3.cumsums)
      seg.list[[paste(model.i, 3)]] <-
        data.frame(chromStart=seg2.chromEnd, chromEnd=last.chromEnd,
                   mean=seg3.means, sample.id=names(profile.list),
                   model.i,
                   row.names=NULL)

      is.seg3 <-
        seg2.chromEnd < small.chromEnd &
          small.chromEnd <= last.chromEnd
      stopifnot(sum(is.seg3) == seg3.bases)
      seg3.mat <- small.bins[is.seg3, , drop=FALSE]
      ##stopifnot(all.equal(colMeans(seg3.mat), seg3.means))

      total.bases <- sum(seg1.bases + seg2.bases + seg3.bases)
      stopifnot(all.equal(total.bases, last.chromEnd - max.chromStart))

      total.loss <- sum(seg1.loss + seg2.loss + seg3.loss)
      model.list[[model.i]] <- data.frame(model.row, total.loss)
    }
    ## Plot the segment means as a reality check.
    seg.df <- do.call(rbind, seg.list)
    ggplot()+
    geom_step(aes(chromStart/1e3, count),
              data=data.frame(norm.df, what="data"),
              color="grey")+
    geom_segment(aes(chromStart/1e3, mean,
                     xend=chromEnd/1e3, yend=mean),
                 data=data.frame(seg.df, what="models"),
                 size=1, color="green")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ model.i, scales="free")

    ## Then plot the peaks only, colored by total cost of the model.
    model.df <- do.call(rbind, model.list)
    model.df$y <- with(model.df, model.i/max(model.i))
    best.model <- model.df[which.min(model.df$total.loss), ]

    ggplot()+
      xlab("position on chromosome (kilobases = kb)")+
      scale_y_continuous("", breaks=NULL)+
    geom_step(aes(chromStart/1e3, count.norm),
              data=data.frame(norm.df, what="data"),
              color="grey")+
    geom_segment(aes(chromStart/1e3, mean.norm,
                     xend=chromEnd/1e3, yend=mean.norm),
                 data=data.frame(bin.df, what="bins"),
                 color="black")+
    geom_segment(aes(peakStart/1e3, 0,
                     xend=peakEnd/1e3, yend=0),
                 data=data.frame(loss.best, what="peak"),
                 color="green")+
    geom_segment(aes(left.chromStart/1e3, y,
                     color=total.loss,
                     xend=right.chromEnd/1e3, yend=y),
                 data=data.frame(model.df, sample.id="peaks"),
                 size=1)+
    geom_text(aes(left.chromStart/1e3, y,
                  label="optimal "),
              data=data.frame(best.model, sample.id="peaks"),
              hjust=1)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ .)

    search.list[[paste(bases.per.bin)]] <-
      data.frame(bases.per.bin, model.df)
    best.list[[paste(bases.per.bin)]] <-
      data.frame(bases.per.bin, best.model)

    peakStart <- best.model$left.chromStart
    peakEnd <- best.model$right.chromEnd
    before.i.list <-
      list(left=best.model$left.cumsum.row-2,
           right=best.model$right.cumsum.row-1)
    for(lr in names(before.i.list)){
      before.i <- before.i.list[[lr]]
      mats <- cumsum.mats[[lr]]
      before.cumsums[[lr]]$count <- mats$count[before.i, ]
    }
  }
  samplefac <- function(L){
    df <- do.call(rbind, L)
    sample.ids <- unique(norm.df$sample.id)
    bases.num <- sort(unique(df$bases.per.bin), decreasing=TRUE)
    levs <- c(paste(sample.ids), paste(bases.num))
    df$sample.id <- factor(df$bases.per.bin, levs)
    df
  }
  search.df <- samplefac(search.list)
  best.df <- samplefac(best.list)
  searchPlot <- 
    ggplot()+
      xlab("position on chromosome (kilobases = kb)")+
      scale_y_continuous("", breaks=NULL)+
    geom_step(aes(chromStart/1e3, count.norm),
              data=data.frame(norm.df, what="data"),
              color="grey")+
    geom_segment(aes(chromStart/1e3, mean.norm,
                     xend=chromEnd/1e3, yend=mean.norm),
                 data=data.frame(bin.df, what="bins"),
                 color="black")+
    geom_vline(aes(xintercept=chromStart/1e3),
               data=optimal,
               color="green")+
    geom_vline(aes(xintercept=chromEnd/1e3),
               data=optimal,
               color="green")+
    geom_segment(aes(left.chromStart/1e3, y,
                     color=total.loss,
                     xend=right.chromEnd/1e3, yend=y),
                 data=search.df,
                 size=1)+
    geom_text(aes(left.chromStart/1e3, y,
                  label="selected"),
              data=best.df,
              hjust=1)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., labeller=function(var, val){
      sub("McGill0", "", val)
    })
  print(searchPlot)
  rbind(optimal=optimal,
        heuristic=data.frame(chromStart=peakStart, chromEnd=peakEnd))

  ## This data set caused a problem in peakStart.
  library(PeakSegDP)
  data(H3K36me3.TDH.other.chunk3.cluster4)
  many <- H3K36me3.TDH.other.chunk3.cluster4
  heuristic.seconds <- system.time({
    heuristic <- multiSampleSegHeuristic(many, 10L)
  })[["elapsed"]]
  min.chromStart <- min(many$chromStart)
  max.chromEnd <- max(many$chromEnd)
  ggplot()+
    geom_step(aes(chromStart/1e3, count), data=many)+
    geom_segment(aes(chromStart/1e3, 0,
                     xend=chromEnd/1e3, yend=0),
                 data=heuristic,
                 color="green",
                 size=1)+
    theme_bw()+
    coord_cartesian(xlim=c(min.chromStart, max.chromEnd)/1e3)+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., scales="free")
})

multiSampleSegOptimal <- function
### Find one peak common to several samples.
(profiles
### List of data.frames with columns chromStart, chromEnd, count, or
### single data.frame with additional column sample.id.
 ){
  if(is.data.frame(profiles)){
    profiles <- split(profiles, profiles$sample.id, drop=TRUE)
  }
  stopifnot(is.list(profiles))
  for(profile.i in seq_along(profiles)){
    df <- profiles[[profile.i]]
    stopifnot(is.data.frame(df))
    stopifnot(is.integer(df$chromStart))
    stopifnot(is.integer(df$chromEnd))
    stopifnot(is.integer(df$count))
    profiles[[profile.i]] <- df[, c("chromStart", "chromEnd", "count")]
  }
  chromStartEnd <-
    .Call("multiSampleSegOptimal_interface",
          profiles,
          PACKAGE="PeakSegDP")
  data.frame(chromStart=chromStartEnd[1],
             chromEnd=chromStartEnd[2])
}

multiSampleSegSome <- structure(function
### Given N samples, find the best overlapping peak in 0, 1, ..., N
### samples.
(profiles,
 bin.factor){
  if(is.data.frame(profiles)){
    profiles <- split(profiles, profiles$sample.id, drop=TRUE)
  }
  stopifnot(is.list(profiles))
  for(profile.i in seq_along(profiles)){
    df <- profiles[[profile.i]]
    stopifnot(is.data.frame(df))
    stopifnot(is.integer(df$chromStart))
    stopifnot(is.integer(df$chromEnd))
    stopifnot(is.integer(df$count))
    profiles[[profile.i]] <- df[, c("chromStart", "chromEnd", "count")]
  }
  seg.models <-
    .Call("multiSampleSegSome_interface",
          profiles,
          PACKAGE="PeakSegDP")
}, ex=function(){
  library(PeakSegJoint)
  data(H3K36me3.TDH.other.chunk1)
  some.counts <-
    subset(H3K36me3.TDH.other.chunk1$counts,
           43100000 < chromEnd & chromStart < 43205000)
  some.regions <- subset(H3K36me3.TDH.other.chunk1$regions,
                         chromStart < 43205000)
  library(ggplot2)
  ann.colors <-
    c(noPeaks="#f6f4bf",
      peakStart="#ffafaf",
      peakEnd="#ff4c4c",
      peaks="#a445ee")
  ggplot()+
    geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                      fill=annotation),
                  alpha=0.5,
                  data=some.regions)+
    scale_fill_manual(values=ann.colors)+
    geom_step(aes(chromStart/1e3, count), data=some.counts)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., scales="free")
  ## Begin R implementation of multiple sample constrained
  ## segmentation heuristic. Input: profiles data.frame.
  profiles <- some.counts
  profile.list <- split(profiles, profiles$sample.id, drop=TRUE)
  max.chromStart <- max(sapply(profile.list, with, chromStart[1]))
  min.chromEnd <- min(sapply(profile.list, with, chromEnd[length(chromEnd)]))
  c(max.chromStart, min.chromEnd)
  bases <- min.chromEnd-max.chromStart
  n.samples <- length(profile.list)
  small.chromEnd <- (max.chromStart+1):min.chromEnd
  small.bins <-
    matrix(NA, bases, n.samples,
           dimnames=list(position=small.chromEnd,
             sample.id=names(profile.list)))
  for(sample.id in names(profile.list)){
    one <- profile.list[[sample.id]]
    bins <- binSum(one, max.chromStart, n.bins=bases)
    stopifnot(bins$chromEnd == small.chromEnd)
    small.bins[, sample.id] <- bins$count
  }
  bin.factor <- 2L
  bases.per.bin <- 1L
  while(bases/bases.per.bin/bin.factor >= 4){
    bases.per.bin <- bases.per.bin * bin.factor
  }
  n.bins <- as.integer(bases %/% bases.per.bin + 1L)
  na.mat <- 
    matrix(NA, n.bins, n.samples,
           dimnames=list(bin=NULL, sample.id=names(profile.list)))
  first.cumsums <- list(count=na.mat)
  bin.list <- list()
  norm.list <- list()
  for(sample.i in seq_along(profile.list)){
    sample.id <- names(profile.list)[sample.i]
    one <- profile.list[[sample.i]]
    max.count <- max(one$count)
    bins <- binSum(one, max.chromStart, bases.per.bin, n.bins)
    stopifnot(n.bins == nrow(bins))
    bins[nrow(bins), "chromEnd"] <- min.chromEnd
    bins$mean <- with(bins, count/(chromEnd-chromStart))
    bins$mean.norm <- bins$mean/max.count
    bin.list[[sample.id]] <- data.frame(sample.id, rbind(bins, NA))
    bases.vec <- with(bins, chromEnd-chromStart)
    stopifnot(bins$count >= 0)
    first.cumsums$count[, sample.i] <- cumsum(bins$count)
    one$count.norm <- one$count/max.count
    norm.list[[sample.i]] <- one
  }
  bin.df <- do.call(rbind, bin.list)
  norm.df <- do.call(rbind, norm.list)
  ## The formula for the optimal Poisson loss 
  ## for 1 segment with d integer data points x_j is
  ## \sum_{j=1}^d m - x_j \log m_j =
  ##   ( \sum_{j=1}^d x_j ) (1-\log m)
  ## where the segment mean m = (\sum x_j)/d,
  OptimalPoissonLoss <- function(mean.value, cumsum.value){
    ifelse(mean.value == 0, 0, cumsum.value * (1-log(mean.value)))
  }
  loss.list <- list()
  peak.list <- list()
  seg.list <- list()
  best.loss.list <- list()
  flat.cumsums <- first.cumsums$count[n.bins, ]
  flat.means <- flat.cumsums/bases
  flat.loss.vec <- OptimalPoissonLoss(flat.means, flat.cumsums)
  best.loss.list[["0"]] <- sum(flat.loss.vec)
  for(seg1.last in 1:(n.bins-2)){
    seg1.cumsums <- first.cumsums$count[seg1.last, ]
    seg1.bases <- seg1.last*bases.per.bin
    seg1.chromEnd <- seg1.bases + max.chromStart
    seg1.means <- seg1.cumsums/seg1.bases
    seg1.loss.vec <- OptimalPoissonLoss(seg1.means, seg1.cumsums)
    seg1.loss <- sum(seg1.loss.vec)
    for(seg2.last in (seg1.last+1):(n.bins-1)){
      cumsum.seg2.end <- first.cumsums$count[seg2.last, ]
      seg2.cumsums <- cumsum.seg2.end-seg1.cumsums
      seg12.bases <- seg2.last*bases.per.bin
      seg2.bases <- seg12.bases - seg1.bases
      seg2.chromEnd <- seg1.chromEnd + seg2.bases
      seg2.means <- seg2.cumsums/seg2.bases
      seg2.loss.vec <- OptimalPoissonLoss(seg2.means, seg2.cumsums)
      seg2.loss <- sum(seg2.loss.vec)
      
      seg3.cumsums <- first.cumsums$count[n.bins, ]-cumsum.seg2.end
      seg3.bases <- bases - seg12.bases
      seg3.means <- seg3.cumsums/seg3.bases

      mean.mat <- rbind(seg1.means, seg2.means, seg3.means)

      this.seg.list <- list()
      for(sample.id in colnames(mean.mat)){
        this.seg.list[[sample.id]] <-
          data.frame(sample.id,
                     chromStart=c(max.chromStart, seg1.chromEnd, seg2.chromEnd,
                                  max.chromStart),
                     chromEnd=c(seg1.chromEnd, seg2.chromEnd, min.chromEnd,
                                min.chromEnd),
                     mean=c(mean.mat[, sample.id], flat.means[[sample.id]]),
                     segments=c(3, 3, 3, 1))
      }
      these.segs <- do.call(rbind, this.seg.list)
      ggplot()+
        scale_color_manual(values=c(data="grey50",
                             bins="black", segments="green"))+
        geom_step(aes(chromStart/1e3, count, color=what),
                  data=data.frame(norm.df, what="data"))+
        geom_segment(aes(chromStart/1e3, mean,
                         xend=chromEnd/1e3, yend=mean,
                         color=what),
                     size=1,
                     data=data.frame(these.segs, what="segments"))+
        geom_segment(aes(chromStart/1e3, mean,
                         xend=chromEnd/1e3, yend=mean,
                         color=what),
                     data=data.frame(bin.df, what="bins"))+
        theme_bw()+
        theme(panel.margin=grid::unit(0, "cm"))+
        facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
          sub("McGill0", "", val)
        })
      
      seg3.loss.vec <- OptimalPoissonLoss(seg3.means, seg3.cumsums)
      seg3.loss <- sum(seg3.loss.vec)
      seg123.loss.vec <- seg1.loss.vec + seg2.loss.vec + seg3.loss.vec

      peak.feasible <- seg1.means < seg2.means & seg2.means > seg3.means
      if(any(peak.feasible)){
        diff.loss.vec <- flat.loss.vec - seg123.loss.vec
        possible.df <- 
          data.frame(flat.loss.vec, seg123.loss.vec,
                     diff.loss.vec, peak.feasible)
        feasible.df <- subset(possible.df, peak.feasible)
        ordered.df <- 
          feasible.df[order(feasible.df$diff.loss.vec, decreasing = TRUE), ]
        for(peaks in 1:nrow(ordered.df)){
          with.peaks <- ordered.df[1:peaks, ]
          with.segs <-
            subset(these.segs,
                   sample.id %in% rownames(with.peaks) & segments==3)
          without.segs <-
            subset(these.segs,
                   !sample.id %in% rownames(with.peaks) & segments==1)
          without.peaks <-
            possible.df[!rownames(possible.df) %in% rownames(with.peaks),]
          with.loss <- with.peaks$seg123.loss.vec
          without.loss <- without.peaks$flat.loss.vec
          total.loss <- sum(with.loss, without.loss)
          loss.list[[paste(peaks)]][[paste(seg1.last, seg2.last)]] <-
            data.frame(seg1.last, seg2.last, peaks, total.loss)
          peak.list[[paste(peaks)]][[paste(seg1.last, seg2.last)]] <- 
            data.frame(sample.id=rownames(with.peaks),
                       chromStart=seg1.last*bases.per.bin+max.chromStart,
                       chromEnd=seg2.last*bases.per.bin+max.chromStart)
          seg.list[[paste(peaks)]][[paste(seg1.last, seg2.last)]] <-
            rbind(without.segs, with.segs)
        }#peaks
      }#any(peak.feasible)
    }#seg2.last
  }#seg1.last
  best.seg.list <- list()
  best.peak.list <- list()
  best.indices.list <- list()
  for(peaks.str in names(loss.list)){
    loss.df <- do.call(rbind, loss.list[[peaks.str]])
    loss.best <- loss.df[which.min(loss.df$total.loss), ]
    best.indices.list[[peaks.str]] <- loss.best
    last.str <- with(loss.best, paste(seg1.last, seg2.last))
    peaks <- as.numeric(peaks.str)
    seg.df <- seg.list[[peaks.str]][[last.str]]

    ggplot()+
      ggtitle(paste0("best model with ", peaks,
                    " peak", ifelse(peaks==1, "", "s")))+
      scale_color_manual(values=c(data="grey50",
                           bins="black", segments="green"))+
      geom_step(aes(chromStart/1e3, count, color=what),
                data=data.frame(norm.df, what="data"))+
      geom_segment(aes(chromStart/1e3, mean,
                       xend=chromEnd/1e3, yend=mean,
                       color=what),
                   size=1,
                   data=data.frame(seg.df, what="segments"))+
      geom_segment(aes(chromStart/1e3, mean,
                       xend=chromEnd/1e3, yend=mean,
                       color=what),
                   data=data.frame(bin.df, what="bins"))+
      theme_bw()+
      theme(panel.margin=grid::unit(0, "cm"))+
      facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
        sub("McGill0", "", val)
      })
    
    best.seg.list[[peaks.str]] <- data.frame(peaks, seg.df)
    best.peak.list[[peaks.str]] <-
      data.frame(peaks, y=peaks*-0.1, peak.list[[peaks.str]][[last.str]])
    best.loss.list[[peaks.str]] <- loss.best$total.loss
  }
  best.peaks <- do.call(rbind, best.peak.list)
  by.sample.loc <-
    split(best.peaks, with(best.peaks, paste(sample.id, chromStart, chromEnd)))
  short.label.list <- list()
  for(sample.loc.name in names(by.sample.loc)){
    sample.loc <- by.sample.loc[[sample.loc.name]]
    peaks.txt <- paste(sample.loc$peaks, collapse=",")
    short.label.list[[sample.loc.name]] <-
      data.frame(sample.loc[1,], peaks.txt)
  }
  short.labels <- do.call(rbind, short.label.list)
  best.peaks$sample.id <- factor(best.peaks$sample.id, names(profile.list))
  sample.counts <- table(best.peaks$sample.id)
  dftype <- function(what, df){
    df$sample.id <- factor(df$sample.id, names(sample.counts))
    data.frame(what, df)
  }
  ggplot()+
    scale_color_manual(values=c(data="grey50",
                         bins="black", peaks="deepskyblue"))+
    geom_step(aes(chromStart/1e3, count.norm, color=what),
              data=dftype("data", norm.df))+
    geom_segment(aes(chromStart/1e3, y,
                     xend=chromEnd/1e3, yend=y,
                     color=what),
                 size=1,
                 data=dftype("peaks", best.peaks))+
    geom_text(aes(chromStart/1e3, y,
                     label=paste0(peaks, " peak",
                       ifelse(peaks==1, "", "s"), " "),
                     color=what),
              hjust=1,
              size=3,
              vjust=0.5,
              data=dftype("peaks", best.peaks))+
    geom_segment(aes(chromStart/1e3, mean.norm,
                     xend=chromEnd/1e3, yend=mean.norm,
                     color=what),
                 data=dftype("bins", bin.df))+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
      sub("McGill0", "", val)
    })
  ggplot()+
    scale_color_manual(values=c(data="grey50",
                         bins="black", peaks="deepskyblue"))+
    geom_step(aes(chromStart/1e3, count, color=what),
              data=dftype("data", norm.df))+
    geom_segment(aes(chromStart/1e3, 0,
                     xend=chromEnd/1e3, yend=0,
                     color=what),
                 size=1,
                 data=dftype("peaks", short.labels))+
    geom_text(aes(chromStart/1e3, 0,
                     label=paste0(peaks.txt, " peaks "),
                     color=what),
                 hjust=1,
              vjust=0,
                 data=dftype("peaks", short.labels))+
    geom_segment(aes(chromStart/1e3, mean,
                     xend=chromEnd/1e3, yend=mean,
                     color=what),
                 data=dftype("bins", bin.df))+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
      sub("McGill0", "", val)
    })
  ## for 0, 1, ..., maxPeaks, run the bin pyramid grid search,
  ## around the peaks found in this first step.
  library(PeakError)
  zoom.peak.list <- list("0"=Peaks())
  zoom.loss.list <-
    list("0"=data.frame(peaks=0, loss=sum(flat.loss.vec)))
  for(peaks.str in names(best.indices.list)){
    loss.best <- best.indices.list[[peaks.str]]
    best.peak.df <- best.peak.list[[peaks.str]]
    samples.with.peaks <- paste(best.peak.df$sample.id)
    sub.norm.df <- subset(norm.df, sample.id %in% samples.with.peaks)
    sub.bin.df <- subset(bin.df, sample.id %in% samples.with.peaks)
    n.samples <- length(samples.with.peaks)
    last.cumsums <- list()
    before.cumsums <- list(left=list(), right=list())
    for(data.type in names(first.cumsums)){
      data.mat <- first.cumsums[[data.type]]
      last.cumsums[[data.type]] <-
        data.mat[nrow(data.mat),][samples.with.peaks]
      before.cumsums$left[[data.type]] <- if(loss.best$seg1.last == 1){
        rep(0, n.samples)
      }else{
        data.mat[loss.best$seg1.last-1,][samples.with.peaks]
      }
      before.cumsums$right[[data.type]] <-
        data.mat[loss.best$seg2.last-1,][samples.with.peaks]
    }
    last.chromEnd <- min.chromEnd
    n.bins.zoom <- bin.factor * 2L
    n.cumsum <- n.bins.zoom + 1L

    ## These will change at the end of each iteration.
    peakStart <- best.peak.df$chromStart[1]
    peakEnd <- best.peak.df$chromEnd[1]
    search.list <- list()
    best.list <- list()
    bases.per.bin.zoom <- bases.per.bin
    while(bases.per.bin.zoom > 1){
      left.chromStart <- peakStart - bases.per.bin.zoom
      right.chromStart <- peakEnd-bases.per.bin.zoom
      bases.per.bin.zoom <- as.integer(bases.per.bin.zoom / bin.factor)

      cumsum.list <- function(chromStart){
        limits <- bases.per.bin.zoom*(0:(n.cumsum-1))+chromStart
        chromStart.vec <- limits[-length(limits)]
        chromEnd.vec <- limits[-1]
        intervals <-
          paste0(chromStart.vec, "-", chromEnd.vec)
        dn <- list(bin=c("before", intervals), sample.id=samples.with.peaks)
        m <- matrix(NA, n.cumsum, n.samples, dimnames=dn)
        list(count=m, 
             chromStart=chromStart.vec, chromEnd=chromEnd.vec)
      }
      cumsum.mats <-
        list(left=cumsum.list(left.chromStart),
             right=cumsum.list(right.chromStart))

      for(sample.i in seq_along(samples.with.peaks)){
        sample.id <- samples.with.peaks[[sample.i]]
        one <- profile.list[[sample.id]]
        lr.list <-
          list(left=binSum(one, left.chromStart,
                 bases.per.bin.zoom, n.bins.zoom),
               right=binSum(one, right.chromStart,
                 bases.per.bin.zoom, n.bins.zoom,
                 empty.as.zero=TRUE))
        for(lr in names(lr.list)){
          lr.bins <- lr.list[[lr]]
          stopifnot(nrow(lr.bins) == n.bins.zoom)
          lr.bases <- with(lr.bins, chromEnd-chromStart)
          lr.before <- before.cumsums[[lr]]
          lr.counts <- list(count=lr.bins$count)
          for(data.type in names(lr.counts)){
            lr.count.vec <-
              c(lr.before[[data.type]][[sample.id]], lr.counts[[data.type]])
            cumsum.mats[[lr]][[data.type]][, sample.id] <-
              cumsum(lr.count.vec)
          }
        }
      }
      possible.grid <- 
        expand.grid(left.cumsum.row=3:n.cumsum, right.cumsum.row=2:n.cumsum)
      possible.grid$left.chromStart <-
        cumsum.mats$left$chromStart[possible.grid$left.cumsum.row-1]
      possible.grid$left.chromEnd <-
        cumsum.mats$left$chromEnd[possible.grid$left.cumsum.row-1]
      possible.grid$right.chromStart <-
        cumsum.mats$right$chromStart[possible.grid$right.cumsum.row-1]
      possible.grid$right.chromEnd <-
        cumsum.mats$right$chromEnd[possible.grid$right.cumsum.row-1]
      feasible.grid <-
        subset(possible.grid,
               left.chromEnd <= right.chromStart &
                 right.chromEnd < last.chromEnd)
      feasible.grid$model.i <- 1:nrow(feasible.grid)
      model.list <- list()
      seg.list <- list()
      sample.loss.list <- list()
      for(model.i in feasible.grid$model.i){
        model.row <- feasible.grid[model.i, ]

        seg1.i <- model.row$left.cumsum.row-1
        seg1.cumsums <- cumsum.mats$left$count[seg1.i, ]
        seg1.chromEnd <- cumsum.mats$left$chromStart[seg1.i]
        seg1.bases <- seg1.chromEnd-max.chromStart
        seg1.means <- seg1.cumsums/seg1.bases
        seg1.loss <- OptimalPoissonLoss(seg1.means, seg1.cumsums)
        seg.list[[paste(model.i, 1)]] <-
          data.frame(chromStart=max.chromStart, chromEnd=seg1.chromEnd,
                     mean=seg1.means, sample.id=samples.with.peaks,
                     model.i,
                     row.names=NULL)

        seg1.mat <-
          small.bins[small.chromEnd <= seg1.chromEnd,
                     samples.with.peaks,
                     drop=FALSE]
        stopifnot(nrow(seg1.mat) == seg1.bases)
        stopifnot(all.equal(as.numeric(colMeans(seg1.mat)),
                            as.numeric(seg1.means)))

        cumsum.seg2.end <-
          cumsum.mats$right$count[model.row$right.cumsum.row, ]
        seg2.cumsums <- cumsum.seg2.end-seg1.cumsums
        seg2.chromEnd <-
          cumsum.mats$right$chromEnd[model.row$right.cumsum.row-1]
        seg2.bases <- seg2.chromEnd-seg1.chromEnd
        seg2.means <- seg2.cumsums/seg2.bases
        seg2.loss <- OptimalPoissonLoss(seg2.means, seg2.cumsums)
        seg.list[[paste(model.i, 2)]] <-
          data.frame(chromStart=seg1.chromEnd, chromEnd=seg2.chromEnd,
                     mean=seg2.means, sample.id=samples.with.peaks,
                     model.i,
                     row.names=NULL)

        is.seg2 <-
          seg1.chromEnd < small.chromEnd &
            small.chromEnd <= seg2.chromEnd
        stopifnot(sum(is.seg2) == seg2.bases)
        seg2.mat <-
          small.bins[is.seg2,
                     samples.with.peaks,
                     drop=FALSE]
        stopifnot(all.equal(as.numeric(colMeans(seg2.mat)),
                            as.numeric(seg2.means)))
        
        seg3.cumsums <- last.cumsums$count-cumsum.seg2.end
        seg3.bases <- last.chromEnd-seg2.chromEnd
        seg3.means <- seg3.cumsums/seg3.bases
        seg3.loss <- OptimalPoissonLoss(seg3.means, seg3.cumsums)
        seg.list[[paste(model.i, 3)]] <-
          data.frame(chromStart=seg2.chromEnd, chromEnd=last.chromEnd,
                     mean=seg3.means, sample.id=samples.with.peaks,
                     model.i,
                     row.names=NULL)

        is.seg3 <-
          seg2.chromEnd < small.chromEnd &
            small.chromEnd <= last.chromEnd
        stopifnot(sum(is.seg3) == seg3.bases)
        seg3.mat <-
          small.bins[is.seg3,
                     samples.with.peaks,
                     drop=FALSE]
        ## stopifnot(all.equal(as.numeric(colMeans(seg3.mat)),
        ##                     as.numeric(seg3.means)))

        total.bases <- sum(seg1.bases + seg2.bases + seg3.bases)
        stopifnot(all.equal(total.bases, last.chromEnd - max.chromStart))

        total.loss <- sum(seg1.loss + seg2.loss + seg3.loss)
        total.loss.vec <- seg1.loss+seg2.loss+seg3.loss
        feasible.vec <- seg1.means < seg2.means & seg2.means > seg3.means
        model.list[[model.i]] <-
          data.frame(model.row, total.loss, feasible=all(feasible.vec))
        sample.loss.list[[model.i]] <-
          data.frame(sample.id=samples.with.peaks, loss=total.loss.vec)
      }
      ## Plot the segment means as a reality check.
      seg.df <- do.call(rbind, seg.list)
      ggplot()+
        geom_step(aes(chromStart/1e3, count),
                  data=data.frame(sub.norm.df, what="data"),
                  color="grey")+
        geom_segment(aes(chromStart/1e3, mean,
                         xend=chromEnd/1e3, yend=mean),
                     data=data.frame(seg.df, what="models"),
                     size=1, color="green")+
        theme_bw()+
        theme(panel.margin=grid::unit(0, "cm"))+
        facet_grid(sample.id ~ model.i, scales="free")

      ## Then plot the peaks only, colored by total cost of the model.
      model.df <- do.call(rbind, model.list)
      model.df$y <- with(model.df, model.i/max(model.i))
      feasible.models <- subset(model.df, feasible)
      best.model <- feasible.models[which.min(feasible.models$total.loss), ]

      ggplot()+
        xlab("position on chromosome (kilobases = kb)")+
        scale_y_continuous("", breaks=NULL)+
        geom_step(aes(chromStart/1e3, count.norm),
                  data=data.frame(sub.norm.df, what="data"),
                  color="grey")+
        geom_segment(aes(chromStart/1e3, mean.norm,
                         xend=chromEnd/1e3, yend=mean.norm),
                     data=data.frame(sub.bin.df, what="bins"),
                     color="black")+
        geom_segment(aes(peakStart/1e3, 0,
                         xend=peakEnd/1e3, yend=0),
                     data=data.frame(loss.best, what="peak"),
                     color="green")+
        scale_linetype_manual(values=c("TRUE"=1, "FALSE"=2))+
        geom_segment(aes(left.chromStart/1e3, y,
                         color=total.loss,
                         linetype=feasible,
                         xend=right.chromEnd/1e3, yend=y),
                     data=data.frame(model.df, sample.id="peaks"),
                     size=1)+
        geom_text(aes(left.chromStart/1e3, y,
                      label="optimal "),
                  data=data.frame(best.model, sample.id="peaks"),
                  hjust=1)+
        theme_bw()+
        theme(panel.margin=grid::unit(0, "cm"))+
        facet_grid(sample.id ~ .)

      search.list[[paste(bases.per.bin.zoom)]] <-
        data.frame(bases.per.bin.zoom, model.df)
      best.list[[paste(bases.per.bin.zoom)]] <-
        data.frame(bases.per.bin.zoom, best.model)

      peakStart <- best.model$left.chromStart
      peakEnd <- best.model$right.chromEnd
      before.i.list <-
        list(left=best.model$left.cumsum.row-2,
             right=best.model$right.cumsum.row-1)
      for(lr in names(before.i.list)){
        before.i <- before.i.list[[lr]]
        mats <- cumsum.mats[[lr]]
        before.cumsums[[lr]]$count <-
          structure(mats$count[before.i,],
                    names=colnames(mats$count))
      }
    }
    samplefac <- function(L){
      df <- do.call(rbind, L)
      sample.ids <- unique(norm.df$sample.id)
      bases.num <- sort(unique(df$bases.per.bin.zoom), decreasing=TRUE)
      levs <- c(paste(sample.ids), paste(bases.num))
      df$sample.id <- factor(df$bases.per.bin.zoom, levs)
      df
    }
    search.df <- samplefac(search.list)
    best.df <- samplefac(best.list)
    searchPlot <- 
      ggplot()+
        xlab("position on chromosome (kilobases = kb)")+
        scale_y_continuous("", breaks=NULL)+
        geom_step(aes(chromStart/1e3, count.norm),
                  data=data.frame(sub.norm.df, what="data"),
                  color="grey")+
        geom_segment(aes(chromStart/1e3, mean.norm,
                         xend=chromEnd/1e3, yend=mean.norm),
                     data=data.frame(sub.bin.df, what="bins"),
                     color="black")+
        geom_segment(aes(left.chromStart/1e3, y,
                         color=total.loss,
                         xend=right.chromEnd/1e3, yend=y),
                     data=search.df,
                     size=1)+
        geom_text(aes(left.chromStart/1e3, y,
                      label="selected"),
                  data=best.df,
                  hjust=1)+
        theme_bw()+
        theme(panel.margin=grid::unit(0, "cm"))+
        facet_grid(sample.id ~ ., labeller=function(var, val){
          sub("McGill0", "", val)
        })
    loss.with.peaks <- sample.loss.list[[best.model$model.i]]
    samples.without.peaks <-
      names(profile.list)[!names(profile.list) %in% samples.with.peaks]
    loss.without.peaks <-
      data.frame(sample.id=samples.without.peaks,
                 loss=flat.loss.vec[samples.without.peaks])
    sample.loss.df <- rbind(loss.with.peaks, loss.without.peaks)
    peaks <- as.numeric(peaks.str)
    zoom.loss.list[[peaks.str]] <- 
      data.frame(peaks, loss=sum(sample.loss.df$loss))
    zoom.peak.list[[peaks.str]] <-
      data.frame(peaks, sample.id=samples.with.peaks,
                 chromStart=peakStart, chromEnd=peakEnd)
  }
  zoom.peaks <- do.call(rbind, zoom.peak.list)
  zoom.loss <- do.call(rbind, zoom.loss.list)
  ggplot()+
    scale_color_manual(values=c(data="grey50",
                         bins="black", peaks="deepskyblue"))+
    geom_step(aes(chromStart/1e3, count.norm, color=what),
              data=dftype("data", norm.df))+
    geom_segment(aes(chromStart/1e3, -peaks*0.1,
                     xend=chromEnd/1e3, yend=-peaks*0.1,
                     color=what),
                 size=1,
                 data=dftype("peaks", zoom.peaks))+
    geom_text(aes(chromStart/1e3, -peaks*0.1,
                  label=paste0(peaks,
                    " peak",
                    ifelse(peaks==1, "", "s"),
                    " "),
                  color=what),
              hjust=1,
              vjust=0.5,
              size=3,
              data=dftype("peaks", zoom.peaks))+
    geom_segment(aes(chromStart/1e3, mean.norm,
                     xend=chromEnd/1e3, yend=mean.norm,
                     color=what),
                 data=dftype("bins", bin.df))+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
      sub("McGill0", "", val)
    })
  ggplot()+
    scale_color_manual(values=c(data="grey50",
                         bins="black", peaks="deepskyblue"))+
    geom_step(aes(chromStart/1e3, count, color=what),
              data=dftype("data", norm.df))+
    geom_segment(aes(chromStart/1e3, 0,
                     xend=chromEnd/1e3, yend=0,
                     color=what),
                 size=1,
                 data=dftype("peaks", zoom.peaks))+
    geom_text(aes(chromStart/1e3, 0,
                     label=paste0(peaks,
                       " peak",
                       ifelse(peaks==1, "", "s"),
                                  " "),
                  color=what),
              hjust=1,
              vjust=0,
              data=dftype("peaks", zoom.peaks))+
    geom_segment(aes(chromStart/1e3, mean,
                     xend=chromEnd/1e3, yend=mean,
                     color=what),
                 data=dftype("bins", bin.df))+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
      sub("McGill0", "", val)
    })
  ggplot()+
    geom_segment(aes(chromStart/1e3, peaks,
                     xend=chromEnd/1e3, yend=peaks),
                 data=zoom.peaks)
  ggplot(zoom.loss, aes(peaks, loss))+
    geom_point()+
    geom_line()
  ## Compute PeakError on this sequence of models.
  regions.by.sample <- split(some.regions, some.regions$sample.id)
  error.list <- list()
  for(peaks.str in names(zoom.peak.list)){
    several.samples <- zoom.peak.list[[peaks.str]]
    peaks.by.sample <- split(several.samples, several.samples$sample.id)
    peaks <- as.numeric(peaks.str)
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
  show.loss <- data.frame(what="Poisson loss", zoom.loss)
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
  rownames(show.error) <- show.error$peaks
  show.loss$errors <- show.error[paste(show.loss$peaks), "errors"]
  show.loss$cummin <- cummin(show.loss$loss)
  compute.loss <- subset(show.loss, loss == cummin)
  exact.df <- with(compute.loss, exactModelSelection(loss, peaks, peaks))
  intercept <- show.loss[as.character(exact.df$peaks), "loss"]
  exact.df$cost <- intercept + exact.df$min.lambda * exact.df$model.complexity
  exact.df$next.cost <- c(exact.df$cost[-1], NA)
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
    geom_abline(aes(slope=peaks, intercept=loss),
                data=show.loss)+
    geom_text(aes(0, loss, label=peaks),
              data=show.loss, hjust=1.5, color="red")+
    ggtitle("model selection: cost = loss_k + lambda*segments_k")
  ## Solve the optimization using grid search.
  L.grid <- with(exact.df,{
    seq(min(max.log.lambda)-1,
        max(min.log.lambda)+1,
        l=100)
  })
  lambda.grid <- exp(L.grid)
  kstar.grid <- sapply(lambda.grid,function(lambda){
    crit <- show.loss$peaks * lambda + show.loss$loss
    picked <- which.min(crit)
    show.loss$peaks[picked]
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
  ## Compute the target interval.
  exact.df$errors <- show.error[paste(exact.df$peaks), "errors"]
  indices <- with(exact.df, {
    largestContinuousMinimum(errors, max.log.lambda-min.log.lambda)
  })
  target.interval <-
    c(exact.df$min.log.lambda[indices$start],
      exact.df$max.log.lambda[indices$end])
  ## Compute the feature matrix for this joint segmentation problem.
  features.by.sample <- list()
  for(sample.id in names(profile.list)){
    ## Compute feature vector for learning using this segmentation
    ## problem.
    sample.counts <- profile.list[[sample.id]]
    too.long <- with(sample.counts, rep(count, chromEnd-chromStart))
    too.long.pos <- with(sample.counts, {
      (chromStart[1]+1):chromEnd[length(chromEnd)]
    })
    stopifnot(length(too.long) == length(too.long.pos))
    keep <- max.chromStart < too.long.pos & too.long.pos <= min.chromEnd
    long <- as.numeric(too.long[keep])
    stopifnot(length(long) == bases)

    feature.vec <-
      c(quartile=quantile(long),
        mean=mean(long),
        sd=sd(long),
        mad=mad(long),
        bases=bases,
        sum=sum(long))
    
    suppressWarnings({
      features.by.sample[[sample.id]] <-
        c(feature.vec,
          `log+1`=log(feature.vec+1),
          log=log(feature.vec),
          log.log=log(log(feature.vec)))
    })
  }
  feature.mat <- do.call(rbind, features.by.sample)
})
