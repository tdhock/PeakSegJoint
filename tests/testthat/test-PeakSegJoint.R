context("PeakSegJoint")

test_that("Step1 C result agrees with R", {
  library(testthat)
  library(PeakSegJoint)
  data(H3K36me3.TDH.other.chunk1)
  lims <- c(43000000, 43200000) # left
  some.sample.ids <- H3K36me3.TDH.other.chunk1$counts$sample.id
  some.counts <-
    subset(H3K36me3.TDH.other.chunk1$counts,
           sample.id %in% some.sample.ids &
           lims[1] < chromEnd & chromStart < lims[2])
  some.regions <- subset(H3K36me3.TDH.other.chunk1$regions,
                         chromStart < lims[2] &
                           sample.id %in% some.sample.ids)
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

  fit <- PeakSegJointHeuristicStep1(some.counts)
  
  ## Begin R implementation of multiple sample constrained
  ## segmentation heuristic. Input: profiles data.frame.
  profiles <- some.counts
  unfilled.profile.list <- split(profiles, profiles$sample.id, drop=TRUE)
  unfilled.chromStart <- max(sapply(unfilled.profile.list, with, chromStart[1]))
  unfilled.chromEnd <-
    min(sapply(unfilled.profile.list, with, chromEnd[length(chromEnd)]))
  unfilled.bases <- unfilled.chromEnd-unfilled.chromStart
  bin.factor <- 2L
  bases.per.bin <- 1L
  while(unfilled.bases/bases.per.bin/bin.factor >= 4){
    bases.per.bin <- bases.per.bin * bin.factor
  }
  n.bins <- as.integer(unfilled.bases %/% bases.per.bin + 1L)

  extra.bases <- n.bins * bases.per.bin - unfilled.bases
  extra.before <- as.integer(extra.bases/2)
  extra.after <- extra.bases - extra.before
  max.chromStart <- unfilled.chromStart-extra.before
  min.chromEnd <- unfilled.chromEnd + extra.after
  profile.list <- list()
  for(sample.id in names(unfilled.profile.list)){
    one.sample <- 
      subset(unfilled.profile.list[[sample.id]],
             unfilled.chromStart < chromEnd &
               chromStart <= unfilled.chromEnd)
    one.sample$chromStart[1] <- unfilled.chromStart
    one.sample$chromEnd[nrow(one.sample)] <- unfilled.chromEnd
    stopifnot(with(one.sample, sum(chromEnd-chromStart)) == unfilled.bases)
    first.row <- last.row <- one.sample[1,]
    first.row$chromStart <- max.chromStart
    first.row$chromEnd <- unfilled.chromStart
    first.row$count <- 0L
    last.row$chromStart <- unfilled.chromEnd
    last.row$chromEnd <- min.chromEnd
    last.row$count <- 0L
    profile.list[[sample.id]] <-
      rbind(first.row, one.sample, last.row)
  }
  bases <- min.chromEnd-max.chromStart
  ## End pre-processing to add zeros.

  expect_identical(fit$seg_start_end[1], max.chromStart)
  expect_identical(fit$seg_start_end[2], min.chromEnd)
  
  ## Small bins are just for testing the computation of the loss
  ## function in the R implementation, and should not be ported to C
  ## code.
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

  ggplot()+
    scale_color_manual(values=c(data="grey50",
                         bins="black", segments="green"))+
    geom_step(aes(chromStart/1e3, count, color=what),
              data=data.frame(norm.df, what="data"))+
    geom_segment(aes(chromStart/1e3, mean,
                     xend=chromEnd/1e3, yend=mean,
                     color=what),
                 data=data.frame(bin.df, what="bins"))+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
      sub("McGill0", "", val)
    })
  
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

  expect_equal(fit$sample_mean_vec, as.numeric(flat.means))
  
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
    model.i <- peaks + 1
    model <- fit$models[[model.i]]
    
    peak.df <- peak.list[[peaks.str]][[last.str]]
    C.start.end <- model$peak_start_end
    R.start.end <- as.integer(peak.df[1, c("chromStart", "chromEnd")])
    expect_equal(C.start.end, R.start.end)

    sample.i <- model$samples_with_peaks_vec + 1
    C.sample.id <- fit$sample.id[sample.i]
    C.sample.sorted <- sort(paste(C.sample.id))
    R.sample.sorted <- sort(paste(peak.df$sample.id))
    expect_identical(C.sample.sorted, R.sample.sorted)

    seg.df <- seg.list[[peaks.str]][[last.str]]
    segs.by.chromStart <- split(seg.df, seg.df$chromStart)
    for(seg.i in seq_along(segs.by.chromStart)){
      seg.i.df <- segs.by.chromStart[[seg.i]]
      rownames(seg.i.df) <- seg.i.df$sample.id
      R.mean.vec <- as.numeric(seg.i.df[C.sample.id, "mean"])
      C.mean.vec <- model[[sprintf("seg%d_mean_vec", seg.i)]]
      expect_equal(C.mean.vec, R.mean.vec)
    }

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
      data.frame(peaks, y=peaks*-0.1, peak.df)
    best.loss.list[[peaks.str]] <- loss.best$total.loss
  }
  R.loss.vec <- as.numeric(best.loss.list)
  C.loss.vec <- sapply(fit$models, "[[", "loss")
  expect_equal(C.loss.vec, R.loss.vec)
  
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
  for(peaks.str in names(best.indices.list)){
    loss.best <- best.indices.list[[peaks.str]]
    best.peak.df <- best.peak.list[[peaks.str]]
    samples.with.peaks <- paste(best.peak.df$sample.id)
    n.samples <- length(samples.with.peaks)
    last.cumsums <- list()
    before.cumsums <- list(left=list(), right=list())
    for(data.type in names(first.cumsums)){
      data.mat <- first.cumsums[[data.type]]
      last.cumsums[[data.type]] <-
        data.mat[nrow(data.mat),][samples.with.peaks]
      before.cumsums$left[[data.type]] <- if(loss.best$seg1.last == 1){
        structure(rep(0, n.samples), names=samples.with.peaks)
      }else{
        data.mat[loss.best$seg1.last-1,][samples.with.peaks]
      }
      before.cumsums$right[[data.type]] <-
        data.mat[loss.best$seg2.last-1,][samples.with.peaks]
    }
    peaks <- as.integer(peaks.str)
    model.i <- peaks+1L
    model <- fit$models[[model.i]]
    C.sample.i <- model$samples_with_peaks_vec + 1L
    C.sample.id <- fit$sample.id[C.sample.i]
    for(lr in names(before.cumsums)){
      R.cumsums <- as.integer(before.cumsums[[lr]]$count[C.sample.id])
      C.cumsums <- model[[paste0(lr, "_cumsum_vec")]]
      expect_equal(C.cumsums, R.cumsums)
    }
  }
  R.last.cumsums <- as.integer(last.cumsums$count[fit$sample.id])
  C.last.cumsums <- fit$last_cumsum_vec
  expect_equal(C.last.cumsums, R.last.cumsums)
  ## for 0, 1, ..., maxPeaks, run the bin pyramid grid search,
  ## around the peaks found in this first step.
})
