ProfileList <- function
### Convert a data.frame or list of profiles to a list that can be
### passed to .Call("PeakSegJointHeuristic...").
(profiles
### List of data.frames with columns chromStart, chromEnd, count, or
### single data.frame with additional column sample.id.
 ){
  if(is.data.frame(profiles)){
    profiles <- as.data.frame(profiles)
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
  profiles
### Named list of data.frames with columns chromStart, chromEnd,
### count.
}

PeakSegJointHeuristicStep1 <- structure(function
### Run the first step of the PeakSegJoint fast heuristic optimization
### algorithm. NB: this function is only for testing the C code
### against the R implementation. For real data see
### PeakSegJointHeuristic.
(profiles,
### List of data.frames with columns chromStart, chromEnd, count, or
### single data.frame with additional column sample.id.
 bin.factor=2L
### Size of bin pyramid. Bigger values result in slower computation.
 ){
  stopifnot(is.numeric(bin.factor))
  stopifnot(length(bin.factor)==1)
  profile.list <- ProfileList(profiles)
  fit <- 
    .Call("PeakSegJointHeuristicStep1_interface",
          profile.list,
          as.integer(bin.factor),
          PACKAGE="PeakSegJoint")
  fit$sample.id <- names(profile.list)
  fit
### List of model fit results, which can be passed to ConvertModelList
### for easier interpretation.
}, ex=function(){
  library(PeakSegJoint)
  library(ggplot2)
  data(H3K36me3.TDH.other.chunk1)
  lims <- c(43000000, 43200000) # left
  some.counts <-
    subset(H3K36me3.TDH.other.chunk1$counts,
           lims[1] < chromEnd & chromStart < lims[2])
  fit <- PeakSegJointHeuristicStep1(some.counts)
  ## Compute bins to show on plot as well.
  profile.list <- split(some.counts, some.counts$sample.id)
  bin.list <- list()
  norm.list <- list()
  for(sample.id in names(profile.list)){
    one <- profile.list[[sample.id]]
    max.count <- max(one$count)
    bins <- binSum(one, fit$bin_start_end[1], fit$bases_per_bin, fit$n_bins)
    stopifnot(fit$n_bins == nrow(bins))
    bins$mean <- with(bins, count/(chromEnd-chromStart))
    bins$mean.norm <- bins$mean/max.count
    stopifnot(bins$count >= 0)
    one$count.norm <- one$count/max.count
    norm.list[[sample.id]] <- one
    bin.list[[sample.id]] <- data.frame(sample.id, bins)
  }
  bin.df <- do.call(rbind, bin.list)
  norm.df <- do.call(rbind, norm.list)
  converted <- ConvertModelList(fit)
  best.peaks <- transform(converted$peaks, y=peaks*-0.1, what="peaks")
  ggplot()+
    scale_color_manual(values=c(data="grey50",
                         peaks="deepskyblue",
                         bins="black", segments="green"))+
    geom_step(aes(chromStart/1e3, count.norm, color=what),
              data=data.frame(norm.df, what="data"))+
    geom_segment(aes(chromStart/1e3, mean.norm,
                     xend=chromEnd/1e3, yend=mean.norm,
                     color=what),
                 data=data.frame(bin.df, what="bins"))+
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
    facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
      sub("McGill0", "", val)
    })
})

PeakSegJointSeveral <- structure(function
### Run PeakSegJoint with several different bin.factor values, keeping
### only the models with lowest Poisson loss for each peak size.
(profiles,
### data.frame or list of them from ProfileList.
 bin.factors=c(2, 3, 5)
### integer vector of suboptimality parameters, bigger values are
### slower.
 ){
  stopifnot(is.numeric(bin.factors))
  stopifnot(length(bin.factors) > 0)
  profile.list <- ProfileList(profiles)
  fit.list <- list()
  for(bf in bin.factors){
    fit.list[[paste(bf)]] <-
      PeakSegJointHeuristic(profile.list, bf)
  }
  best.fit <- fit.list[[1]]
  for(model.i in seq_along(best.fit$models)){
    for(fit in fit.list[-1]){
      model <- fit$models[[model.i]]
      if(model$loss < best.fit$models[[model.i]]$loss){
        best.fit$models[[model.i]] <- model
      }
    }
  }
  best.fit
### Model fit list as is returned by PeakSegJointHeuristic.
}, ex=function(){
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
})

PeakSegJointHeuristic <- structure(function
### Run the PeakSegJoint fast heuristic optimization algorithm, which
### gives an approximate solution to a multi-sample Poisson maximum
### likelihood segmentation problem. Given S samples, this function
### computes a sequence of S+1 PeakSegJoint models, with 0, ..., S
### samples with an overlapping peak (maximum of one peak per sample).
(profiles,
### List of data.frames with columns chromStart, chromEnd, count, or
### single data.frame with additional column sample.id.
 bin.factor=2L
### Size of bin pyramid. Bigger values result in slower computation.
 ){
  stopifnot(is.numeric(bin.factor))
  stopifnot(length(bin.factor)==1)
  profile.list <- ProfileList(profiles)
  fit <- 
    .Call("PeakSegJointHeuristic_interface",
          profile.list,
          as.integer(bin.factor),
          PACKAGE="PeakSegJoint")
  fit$sample.id <- names(profile.list)
  fit
### List of model fit results, which can be passed to ConvertModelList
### for easier interpretation.
}, ex=function(){
  library(PeakSegJoint)
  data(H3K36me3.TDH.other.chunk1)
  lims <- c(43000000, 43200000) # left
  some.counts <-
    subset(H3K36me3.TDH.other.chunk1$counts,
           lims[1] < chromEnd & chromStart < lims[2])
  library(microbenchmark)
  microbenchmark(fit={
    fit <- PeakSegJointHeuristic(some.counts)
  }, fit.convert={
    fit <- PeakSegJointHeuristic(some.counts)
    converted <- ConvertModelList(fit)
  }, times=10)
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
  library(ggplot2)
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
    facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
      sub("McGill0", "", val)
    })
  ggplot()+
    geom_segment(aes(chromStart/1e3, peaks,
                     xend=chromEnd/1e3, yend=peaks),
                 data=best.peaks)
  ggplot(converted$loss, aes(peaks, loss))+
    geom_point()+
    geom_line()
  ## PeakSegJointHeuristic can also be used as a fast approximate
  ## solver for the PeakSeg model with just 1 peak.
  some.counts <-
    subset(H3K36me3.TDH.other.chunk1$counts,
           43379893 < chromEnd)
  profile.list <- split(some.counts, some.counts$sample.id)
  require(PeakSegDP)
  sample.peak.list <- list()
  for(sample.id in names(profile.list)){
    sample.counts <- profile.list[[sample.id]]
    fit <- PeakSegJointHeuristic(sample.counts)
    heuristic <- fit$models[[2]]$peak_start_end
    dp <- PeakSegDP(sample.counts, maxPeaks=1L)
    dp.peak <- dp$peaks[["1"]]
    sample.peak.list[[sample.id]] <-
      data.frame(sample.id,
                 y=-0.1 * (1:2) * max(sample.counts$count),
                 algorithm=c("heuristic", "cDPA"),
                 chromStart=c(heuristic[1], dp.peak$chromStart),
                 chromEnd=c(heuristic[2], dp.peak$chromEnd))
  }
  sample.peaks <- do.call(rbind, sample.peak.list)
  ggplot()+
    scale_color_discrete(limits=c("heuristic", "cDPA"))+
    geom_step(aes(chromStart/1e3, count),
              data=some.counts,
              color="grey50")+
    geom_segment(aes(chromStart/1e3, y,
                     color=algorithm,
                     xend=chromEnd/1e3, yend=y),
                 data=sample.peaks,
                 size=2)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., scales="free", labeller=function(var, val){
      sub("McGill0", "", val)
    })
})

ConvertModelList <- function
### Convert a model list from the non-repetitive format that we get
### from the C code to the repetitive format that is more useful for
### plotting.
(model.list
### Value of PeakSegJointHeuristicStep1(...).
 ){
  seg1.chromStart <- model.list$data_start_end[1]
  seg3.chromEnd <- model.list$data_start_end[2]
  seg.list <- list()
  loss.list <- list()
  peak.list <- list()
  for(model.i in seq_along(model.list$models)){
    model <- model.list$models[[model.i]]
    peaks <- model.i-1
    peaks.str <- paste(peaks)
    loss <- model$loss
    loss.list[[peaks.str]] <- data.frame(peaks, loss)
    if(is.finite(loss)){
      sample.i.vec <- model$samples_with_peaks_vec + 1
      samples.with.peaks <- model.list$sample.id[sample.i.vec]
      samples.without.peaks <- model.list$sample.id[-sample.i.vec]
      peakStart <- model$peak_start_end[1]
      peakEnd <- model$peak_start_end[2]
      if(peaks > 0){
        seg.list[[paste(peaks, 1)]] <-
          data.frame(peaks,
                     sample.id=samples.with.peaks,
                     chromStart=seg1.chromStart,
                     chromEnd=peakStart,
                     mean=model$seg1_mean_vec)
        peak.list[[peaks.str]] <- seg.list[[paste(peaks, 2)]] <- 
          data.frame(peaks,
                     sample.id=samples.with.peaks,
                     chromStart=peakStart,
                     chromEnd=peakEnd,
                     mean=model$seg2_mean_vec)
        seg.list[[paste(peaks, 3)]] <- 
          data.frame(peaks,
                     sample.id=samples.with.peaks,
                     chromStart=peakEnd,
                     chromEnd=seg3.chromEnd,
                     mean=model$seg3_mean_vec)
      }
      if(length(samples.without.peaks)){
        seg.list[[paste(peaks, "flat")]] <- 
          data.frame(peaks,
                     sample.id=samples.without.peaks,
                     chromStart=seg1.chromStart,
                     chromEnd=seg3.chromEnd,
                     mean=model.list$sample_mean_vec[-sample.i.vec])
      }
    }
  }
  list(peaks=do.call(rbind, peak.list),
       segments=do.call(rbind, seg.list),
       loss=do.call(rbind, loss.list))
### List of peaks, segments, loss.
}  
