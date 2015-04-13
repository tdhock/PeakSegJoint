ProfileList <- function(profiles){
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
}

PeakSegJointHeuristicStep1 <- structure(function
### Run the first step of the PeakSegJoint fast heuristic optimization
### algorithm, for testing the results of the C code against the R
### implementation.
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
    bins <- binSum(one, fit$seg_start_end[1], fit$bases_per_bin, fit$n_bins)
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

PeakSegJointHeuristic <- structure(function
### Run the PeakSegJoint fast heuristic optimization algorithm, which
### attempts to find the 1 or 3 segment models with maximum Poisson
### likelihood, across several profiles.
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
  })
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
})

ConvertModelList <- function
### Convert a model list from the non-repetitive format that we get
### from the C code to the repetitive format that is more useful for
### plotting.
(model.list
### Value of PeakSegJointHeuristicStep1(...).
 ){
  seg1.chromStart <- model.list$seg_start_end[1]
  seg3.chromEnd <- model.list$seg_start_end[2]
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
