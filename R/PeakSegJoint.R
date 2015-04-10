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

PeakSegJointHeuristicStep1 <- function
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
}

ConvertModelList <- function(model.list){
  seg1.chromStart <- model.list$seg_start_end[1]
  seg3.chromEnd <- model.list$seg_start_end[2]
  seg.list <- list()
  loss.list <- list()
  peak.list <- list()
  cumsum.list <- list()
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
        cumsum.list[[peaks.str]] <-
          data.frame(peaks,
                     sample.id=samples.with.peaks,
                     left=model$left_cumsum_vec,
                     right=model$right_cumsum_vec,
                     last=model$last_cumsum_vec)
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
       cumsum=do.call(rbind, cumsum.list),
       loss=do.call(rbind, loss.list))
}  
