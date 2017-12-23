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
    sample.path <- if("sample.group" %in% names(profiles)){
      with(profiles, paste0(sample.group, "/", sample.id))
    }else{
      profiles$sample.id
    }
    profiles <- split(profiles, sample.path, drop=TRUE)
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
### algorithm. This is the GridSearch subroutine of the JointZoom algorithm
### in arXiv:1506.01286. NB: this function is only for testing the C code
### against the R implementation (search tests/testthat/*.R for Step1).
### For real data see
### PeakSegJointSeveral.
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
  data(H3K36me3.TDH.other.chunk1, envir=environment())
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

  if(require(ggplot2) && interactive()){

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
      facet_grid(sample.id ~ ., scales="free")

  }
  
})

PeakSegJointHeuristicStep2 <- function
### Run the first and second steps of the PeakSegJoint fast heuristic
### optimization algorithm. Step2 the SearchNearPeak subroutine
### described in the JointZoom Algorithm of arXiv:1506.01286, and it
### is guaranteed to return feasible segmentations (seg1 < seg2 >
### seg3). NB: this function is only for testing the C code against
### the R implementation (search tests/testthat/*.R files for
### Step2). For real data see PeakSegJointSeveral.
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
    .Call("PeakSegJointHeuristicStep2_interface",
          profile.list,
          as.integer(bin.factor),
          PACKAGE="PeakSegJoint")
  fit$sample.id <- names(profile.list)
  fit
### List of model fit results, which can be passed to ConvertModelList
### for easier interpretation.
}

### Run the PeakSegJoint heuristic segmentation algorithm with several
### different bin.factor values, keeping only the models with lowest
### Poisson loss for each peak size. This algorithm gives an
### approximate solution to the following multi-sample constrained
### maximum likelihood segmentation problem. If there are S samples
### total, we look for the most likely common peak in \eqn{s\in{0, ..., S}}
### samples. We solve the equivalent minimization problem using the
### Poisson loss seg.mean - count.data * log(seg.mean), from the first
### base to the last base of profiles. The optimization variables are
### the segment means, of which there can be either 1 value (no peak)
### or 3 values (peak) in each sample. If there are 3 segments then
### two constraints are applied: (1) the changes in mean must occur at
### the same position in each sample, and (2) the changes must be up
### and then down (mean1 < mean2 > mean3).
PeakSegJointSeveral <- structure(function
(profiles,
### data.frame or list of them from ProfileList.
 bin.factors=2:7
### integer vector of optimization parameters >= 2. Larger values are
### slower. Using more values is slower since it tells the algorithm
### to search more of the model space, yielding solution which is
### closer to the global optimum.
 ){
  stopifnot(is.numeric(bin.factors))
  stopifnot(length(bin.factors) > 0)
  profile.list <- ProfileList(profiles)
  fit.list <- list()
  for(bf in bin.factors){
    fit.list[[paste(bf)]] <- tryCatch({
      PeakSegJointHeuristic(profile.list, bf)
    }, error=function(e){
      NULL
    })
  }
  if(length(fit.list) == 0){
    stop("No computable models")
  }
  best.fit <- fit.list[[1]]
  for(fit in fit.list[-1]){
    fit.loss <- sapply(fit$models, "[[", "loss")
    best.loss <- sapply(best.fit$models, "[[", "loss")
    to.copy <- fit.loss < best.loss
    best.fit$models[to.copy] <- fit$models[to.copy]
  }
  best.fit
### List of model fit results, which can be passed to ConvertModelList
### for easier interpretation.
}, ex=function(){

  library(PeakSegJoint)
  data(H3K4me3.TDH.other.chunk8, envir=environment())
  bf.vec <- c(2, 3, 5)
  fit.list <-
    list(several=PeakSegJointSeveral(H3K4me3.TDH.other.chunk8, bf.vec))
  for(bf in bf.vec){
    fit.list[[paste(bf)]] <-
      PeakSegJointHeuristicStep2(H3K4me3.TDH.other.chunk8, bf)
  }
  loss.list <- list()
  segs.by.peaks.fit <- list()
  for(fit.name in names(fit.list)){
    fit <- fit.list[[fit.name]]
    loss.list[[fit.name]] <- sapply(fit$models, "[[", "loss")
    converted <- ConvertModelList(fit)
    segs.by.peaks <- with(converted, split(segments, segments$peaks))
    for(peaks in names(segs.by.peaks)){
      model.segs <- segs.by.peaks[[peaks]]
      if(is.data.frame(model.segs)){
        segs.by.peaks.fit[[peaks]][[fit.name]] <-
          data.frame(fit.name, model.segs)
      }
    }
  }
  do.call(rbind, loss.list)

  segs1 <- do.call(rbind, segs.by.peaks.fit[["10"]])
  breaks1 <- subset(segs1, min(chromStart) < chromStart)
  if(interactive() && require(ggplot2)){
    ggplot()+
      ggtitle(paste("PeakSegJointSeveral runs PeakSegJointHeuristic",
                    "and keeps only the most likely model"))+
      geom_step(aes(chromStart/1e3, count),
                color="grey50",
                data=H3K4me3.TDH.other.chunk8)+
      geom_vline(aes(xintercept=chromStart/1e3),
                 data=breaks1,
                 color="green",
                 linetype="dashed")+
      geom_segment(aes(chromStart/1e3, mean,
                       xend=chromEnd/1e3, yend=mean),
                   size=1,
                   color="green",
                   data=segs1)+
      theme_bw()+
      theme(panel.margin=grid::unit(0, "cm"))+
      facet_grid(sample.id ~ fit.name, scales="free")
  }

  segs.by.peaks <- list()
  for(peaks in 8:10){
    segs.by.peaks[[paste(peaks)]] <-
      data.frame(peaks, segs.by.peaks.fit[[paste(peaks)]][["several"]])
  }
  segs <- do.call(rbind, segs.by.peaks)
  breaks <- subset(segs, min(chromStart) < chromStart)
  if(interactive() && require(ggplot2)){
    ggplot()+
      ggtitle("PeakSegJoint models with 8-10 peaks")+
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
      facet_grid(sample.id ~ peaks, scales="free")
  }
  
})

### Run the PeakSegJoint fast heuristic optimization algorithm, which
### gives an approximate solution to a multi-sample Poisson maximum
### likelihood segmentation problem. Given S samples, this function
### computes a sequence of S+1 PeakSegJoint models, with 0, ..., S
### samples with an overlapping peak (maximum of one peak per
### sample). This solver runs steps 1-3, and Step3 checks if there are
### any more likely models in samples with peak locations which are
### the same as all the models detected in Step2. This is guaranteed
### as of 24 July 2015 to return a feasible segmentation (seg1 < seg2
### > seg3). NB: this function is mostly for internal testing purposes
### (search tests/testthat/*.R for 'Heuristic('). For real data use
### PeakSegJointSeveral.
PeakSegJointHeuristic <- structure(function
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
  data(H3K36me3.TDH.other.chunk1, envir=environment())
  lims <- c(43000000, 43200000) # left
  some.counts <-
    subset(H3K36me3.TDH.other.chunk1$counts,
           lims[1] < chromEnd & chromStart < lims[2])
  fit <- PeakSegJointHeuristic(some.counts)
  converted <- ConvertModelList(fit)
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
  
  if(interactive() && require(ggplot2)){
    
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
      facet_grid(sample.id ~ ., scales="free")
  
    ggplot(converted$loss, aes(peaks, loss))+
      geom_point()+
      geom_line()

  }

})

### Run the PeakSegJointFaster heuristic optimization algorithm, which
### gives an approximate solution to a multi-sample Poisson maximum
### likelihood segmentation problem. Given S samples, this function
### computes a sequence of S+1 PeakSegJoint models, with 0, ..., S
### samples with an overlapping peak (maximum of one peak per
### sample). 
PeakSegJointFasterOne <- structure(function
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
    .Call("PeakSegJointFaster_interface",
          profile.list,
          as.integer(bin.factor),
          PACKAGE="PeakSegJoint")
  fit$sample.id <- names(profile.list)
  fit
### List of model fit results, see examples to see how to use it.
}, ex=function(){

  library(PeakSegJoint)
  data(H3K36me3.TDH.other.chunk1, envir=environment())

  some.counts <- subset(
    H3K36me3.TDH.other.chunk1$counts,
    43000000 < chromEnd &
    chromStart < 43200000 &
    sample.id %in% c("McGill0023", "McGill0022", "McGill0016", "McGill0013"))

  id.df <- unique(some.counts[, c("cell.type", "sample.id")])
  group.list <- split(paste(id.df$sample.id), id.df$cell.type, drop=TRUE)

  loss.df.list <- list()
  fit.list <- list()
  for(bin.factor in 2:7){
    fit.fast <- PeakSegJointFasterOne(some.counts, bin.factor)
    fit.fast$min.loss <- sum(fit.fast$peak_loss_vec)
    fit.fast$sample.loss.diff.vec <- sort(with(fit.fast, structure(
      peak_loss_vec-flat_loss_vec, names=sample.id)))
    fit.fast$group.loss.diff.vec <- sort(sapply(group.list, function(sid.vec){
      sum(fit.fast$sample.loss.diff.vec[sid.vec])
    }))
    fit.fast$sample.loss.vec <- with(fit.fast, structure(
      sum(flat_loss_vec)+cumsum(c(0, sample.loss.diff.vec)),
      names=paste0(0:length(sample.loss.diff.vec), "samples")))
    fit.fast$group.loss.vec <- with(fit.fast, structure(
      sum(flat_loss_vec)+cumsum(c(0, group.loss.diff.vec)),
      names=paste0(0:length(group.loss.diff.vec), "groups")))
    loss.df.list[[paste(bin.factor)]] <- with(fit.fast, data.frame(
      bin.factor,
      loss=sample.loss.vec,
      peaks=0:length(sample.loss.diff.vec)))
    fit.list[[paste(bin.factor)]] <- fit.fast
  }
  loss.df <- do.call(rbind, loss.df.list)
  fit.best <- fit.list[[which.min(sapply(fit.list, "[[", "min.loss"))]]

  norm.list <- list()
  profile.list <- split(some.counts, some.counts$sample.id, drop=TRUE)
  for(sample.id in names(profile.list)){
    one <- profile.list[[sample.id]]
    max.count <- max(one$count)
    one$count.norm <- one$count/max.count
    norm.list[[sample.id]] <- one
  }
  norm.df <- do.call(rbind, norm.list)

  if(interactive() && require(ggplot2)){
    
    peaks.df.list <- list()
    for(n.samples in 1:length(fit.best$sample.loss.diff.vec)){
      peaks.df.list[[paste(n.samples)]] <- with(fit.best, data.frame(
        samples=n.samples,
        sample.id=names(sample.loss.diff.vec)[1:n.samples],
        chromStart=peak_start_end[1],
        chromEnd=peak_start_end[2]))
    }
    peaks <- do.call(rbind, peaks.df.list)
    best.peaks <- transform(peaks, y=samples*-0.1, what="peaks")
    ggplot()+
      ggtitle("model for each sample")+
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
                    label=paste0(samples, " sample",
                      ifelse(samples==1, "", "s"), " "),
                    color=what),
                hjust=1,
                size=3,
                vjust=0.5,
                data=best.peaks)+
      theme_bw()+
      theme(panel.margin=grid::unit(0, "cm"))+
      facet_grid(sample.id ~ ., scales="free")

    ## same thing but for each group.
    peaks.df.list <- list()
    for(n.groups in 1:length(fit.best$group.loss.diff.vec)){
      group.vec <- names(fit.best$group.loss.diff.vec[1:n.groups])
      meta.df <- do.call(rbind, lapply(group.vec, function(cell.type){
        data.frame(cell.type, sample.id=group.list[[cell.type]])
      }))
      peaks.df.list[[paste(n.groups)]] <- with(fit.best, data.frame(
        groups=n.groups,
        meta.df,
        chromStart=peak_start_end[1],
        chromEnd=peak_start_end[2]))
    }
    peaks <- do.call(rbind, peaks.df.list)
    best.peaks <- transform(peaks, y=groups*-0.1, what="peaks")
    ggplot()+
      ggtitle("model for each group")+
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
                    label=paste0(groups, " group",
                      ifelse(groups==1, "", "s"), " "),
                    color=what),
                hjust=1,
                size=3,
                vjust=0.5,
                data=best.peaks)+
      theme_bw()+
      theme(panel.margin=grid::unit(0, "cm"))+
      facet_grid(sample.id + cell.type ~ ., scales="free")
  
    min.df <- subset(loss.df, peaks==max(peaks))
    ggplot()+
      geom_line(aes(peaks, loss, group=bin.factor), data=loss.df)+
      geom_text(aes(peaks, loss, label=bin.factor), data=min.df, hjust=0)

    if(require(microbenchmark)){

      N.samples.vec <- 10^seq(1, 3, by=0.5)
      max.N <- max(N.samples.vec)
      N.bases <- 10
      rmat <- function(Nr, Nc, mu){
        matrix(rpois(Nr*Nc, mu), Nr, Nc)
      }
      set.seed(1)
      big.mat <- cbind(
        rmat(max.N, N.bases, 5),
        rmat(max.N, N.bases, 10),
        rmat(max.N, N.bases, 5))
      big.df <- data.frame(
        sample.id=as.integer(row(big.mat)),
        chromStart=as.integer(col(big.mat)-1),
        chromEnd=as.integer(col(big.mat)),
        count=as.integer(big.mat))
      full.list <- ProfileList(big.df)
      time.df.list <- list()
      for(N.samples in N.samples.vec){
        partial.list <- full.list[1:N.samples]
        result <- microbenchmark(
          Heuristic=PeakSegJointHeuristic(partial.list, 2L),
          Faster=PeakSegJointFasterOne(partial.list, 2L),
          times=2L)
        time.df.list[[paste(N.samples)]] <- data.frame(
          N.samples,
          result)
      }
      time.df <- do.call(rbind, time.df.list)

      ggplot()+
        geom_point(aes(
          N.samples, time/1e9, color=expr),
          data=time.df)+
        scale_x_log10()+
        scale_y_log10("seconds")
      
    }

  }

})

### Run the PeakSegJointFaster heuristic optimization algorithm, for
### several bin.factor parameter values, keeping only the most likely
### model found. This gives an approximate solution to a multi-sample
### Poisson maximum likelihood segmentation problem. Given S samples,
### this function computes a sequence of S+1 PeakSegJoint models, with
### 0, ..., S samples with an overlapping peak (maximum of one peak
### per sample). It also computes for G groups, the seq of G+1
### models, with 0, ..., G groups with an overlapping peak.
PeakSegJointFaster <- structure(function
(profiles,
### data.frame with columns sample.id, sample.group, chromStart,
### chromEnd, count.
 bin.factor=2:7
### Size of bin pyramid. Bigger values result in slower computation.
 ){
  stopifnot(is.numeric(bin.factor))
  stopifnot(length(bin.factor)>0)
  profile.list <- ProfileList(profiles)
  bad.names <- grep("/", names(profile.list), invert=TRUE, value=TRUE)
  if(length(bad.names)){
    print(bad.names)
    stop("no sample.group column in profile data but this is needed for PeakSegJointFaster")
  }
  sample.group <- sub("/.*", "", names(profile.list))
  sample.id <- sub(".*/", "", names(profile.list))
  id.df <- data.frame(
    sample.group, sample.id,
    sample.path=paste0(sample.group, "/", sample.id))
  group.list <- split(paste(id.df$sample.path), id.df$sample.group, drop=TRUE)
  fit.list <- list()
  for(bin.factor in 2:7){
    fit.fast <- PeakSegJointFasterOne(profile.list, bin.factor)
    fit.fast$min.loss <- sum(fit.fast$peak_loss_vec)
    fit.list[[paste(bin.factor)]] <- fit.fast
  }
  fit.best <- fit.list[[which.min(sapply(fit.list, "[[", "min.loss"))]]
  rownames(fit.best$mean_mat) <- fit.best$sample.id
  is.feasible.vec <- with(fit.best, {
    mean_mat[,1] < mean_mat[,2] & mean_mat[,2] > mean_mat[,3]
  })
  feasible.name.vec <- names(is.feasible.vec)[is.feasible.vec]
  fit.best$group.list <- group.list
  fit.best$sample.loss.diff.vec <- sort(with(fit.best, structure(
    peak_loss_vec-flat_loss_vec, names=sample.id))[feasible.name.vec])
  fit.best$group.loss.diff.vec <- sort(sapply(group.list, function(sid.vec){
    sum(fit.best$sample.loss.diff.vec[sid.vec])
  }))
  fit.best$sample.loss.vec <- with(fit.best, structure(
    sum(flat_loss_vec)+cumsum(c(0, sample.loss.diff.vec)),
    names=paste0(0:length(sample.loss.diff.vec), "samples")))
  fit.best$group.loss.vec <- with(fit.best, structure(
    sum(flat_loss_vec)+cumsum(c(0, group.loss.diff.vec)),
    names=paste0(0:length(group.loss.diff.vec), "groups")))
  fit.best$group.modelSelection <- with(fit.best, {
    penaltyLearning::modelSelection(data.frame(
      loss=group.loss.vec,
      complexity=0:length(group.loss.diff.vec)))
  })
  fit.best$sample.modelSelection <- with(fit.best, {
    penaltyLearning::modelSelection(data.frame(
      loss=sample.loss.vec,
      complexity=0:length(sample.loss.diff.vec)))
  })
  fit.best
### List of model fit results.
}, ex=function(){

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

})

ConvertModelList <- function
### Convert a model list from the non-repetitive format that we get
### from the C code to the repetitive format that is more useful for
### plotting.
(model.list
### List from PeakSegJointHeuristic(...) or PeakSegJointSeveral(...).
 ){
  seg1.chromStart <- model.list$data_start_end[1]
  seg3.chromEnd <- model.list$data_start_end[2]
  bases <- seg3.chromEnd - seg1.chromStart
  seg.list <- list()
  loss.list <- list()
  peak.list <- list()
  n.samples <- length(model.list$sample.id)
  for(model.i in seq_along(model.list$models)){
    model <- model.list$models[[model.i]]
    peaks <- model.i-1
    peaks.str <- paste(peaks)
    loss <- model$loss
    segments.vec <- rep(1, n.samples)
    if(0 < peaks){
      segments.vec[1:peaks] <- 3
    }
    under.sqrt <- 1.1 + log(bases/segments.vec)
    in.square <- 1 + 4 * sqrt(under.sqrt)
    model.complexity <- sum(segments.vec * in.square * in.square)
    loss.list[[peaks.str]] <- data.frame(peaks, loss, model.complexity)
    if(is.finite(loss)){
      sample.i.vec <- model$samples_with_peaks_vec + 1
      has.peak <- seq_along(model.list$sample.id) %in% sample.i.vec
      samples.with.peaks <- model.list$sample.id[sample.i.vec]
      samples.without.peaks <- model.list$sample.id[!has.peak]
      parsePath <- function(sample.path){
        sample.id <- sub(".*/", "", sample.path)
        maybe.group <- sub("/.*", "", sample.path)
        sample.group <- ifelse(maybe.group == sample.path, NA, maybe.group)
        if(all(is.na(sample.group))){
          data.frame(peaks, sample.id)
        }else{
          data.frame(peaks, sample.id, sample.group)
        }
      }
      peakStart <- model$peak_start_end[1]
      peakEnd <- model$peak_start_end[2]
      if(peaks > 0){
        meta <- parsePath(samples.with.peaks)
        seg.list[[paste(peaks, 1)]] <-
          data.frame(meta,
                     chromStart=seg1.chromStart,
                     chromEnd=peakStart,
                     mean=model$seg1_mean_vec)
        peak.list[[peaks.str]] <- seg.list[[paste(peaks, 2)]] <- 
          data.frame(meta,
                     chromStart=peakStart,
                     chromEnd=peakEnd,
                     mean=model$seg2_mean_vec)
        seg.list[[paste(peaks, 3)]] <- 
          data.frame(meta,
                     chromStart=peakEnd,
                     chromEnd=seg3.chromEnd,
                     mean=model$seg3_mean_vec)
      }
      if(length(samples.without.peaks)){
        seg.list[[paste(peaks, "flat")]] <- 
          data.frame(parsePath(samples.without.peaks),
                     chromStart=seg1.chromStart,
                     chromEnd=seg3.chromEnd,
                     mean=model.list$sample_mean_vec[!has.peak])
      }
    }
  }
  info <- 
    list(peaks=do.call(rbind, peak.list),
         segments=do.call(rbind, seg.list),
         loss=do.call(rbind, loss.list))
  seg.size <- with(info$segments, chromEnd - chromStart)
  too.small <- seg.size < 1
  if(any(too.small)){
    print(info$segments[too.small, ])
    stop("solver reported segment less than 1 base")
  }
  info$loss$cummin <- cummin(info$loss$loss)
  some.loss <- subset(info$loss, loss == cummin)
  ## what to do when models have the same loss but different number of
  ## peaks? keep the model with min peaks:
  cummin.reduced <- c(TRUE, diff(some.loss$cummin) < 0)
  ## keep the model with max peaks:
  cummin.reduced <- rev(c(TRUE, diff(rev(some.loss$cummin)) > 0))
  decreasing.loss <- some.loss[cummin.reduced, ]
  info$modelSelection <- penaltyLearning::modelSelection(
    decreasing.loss, "loss", "peaks")
  info
### List of data.frames: segments has 1 row for each segment mean,
### sample, and model size (peaks, sample.id, sample.group,
### chromStart, chromEnd, mean); peaks is the same kind of data.frame
### as segments, but with only the second/peak segments; loss has one
### row for each model size; modelSelection has one row for each model
### size that can be selected, see exactModelSelection.
}  
