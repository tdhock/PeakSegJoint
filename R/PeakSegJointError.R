PeakSegJointError <- structure(function
### Compute number of incorrect regions for every PeakSegJoint model.
(converted,
### Result of ConvertModelList.
 problem.regions
### data.frame of annotated region labels.
 ){
  getID <- function(df){
    if("sample.group" %in% names(df)){
      with(df, paste0(sample.group, "/", sample.id))
    }else{
      paste(df$sample.id)
    }
  }
  regions.by.sample <- split(problem.regions, getID(problem.regions))
  peaks.by.peaks <- with(converted, if(is.null(peaks)){
    list()
  }else{
    split(peaks, peaks$peaks)
  })
  error.by.peaks <- list()
  error.regions.list <- list()
  for(peaks.str in c("0", names(peaks.by.peaks))){
    peaks.by.sample <- if(peaks.str=="0"){
      list()
    }else{
      model.peaks <- peaks.by.peaks[[peaks.str]]
      split(model.peaks, getID(model.peaks))
    }
    error.by.sample <- list()
    for(sample.path in names(regions.by.sample)){
      sample.regions <- regions.by.sample[[sample.path]]
      sample.peaks <- if(sample.path %in% names(peaks.by.sample)){
        peaks.by.sample[[sample.path]]
      }else{
        Peaks()
      }
      model.error <- PeakErrorChrom(sample.peaks, sample.regions)
      sample.id <- sub(".*/", "", sample.path)
      sample.group <- sub("/.*", "", sample.path)
      error.by.sample[[sample.path]] <-
        data.frame(sample.id, sample.group,
                   peaks=as.integer(peaks.str),
                   model.error,
                   row.names=NULL)
    }
    model.error <- do.call(rbind, error.by.sample)
    rownames(model.error) <- NULL
    error.regions.list[[peaks.str]] <- model.error
    error.by.peaks[[peaks.str]] <- with(model.error, {
      data.frame(peaks=as.integer(peaks.str),
                 fp=sum(fp),
                 fn=sum(fn),
                 tp=sum(tp),
                 possible.fp=sum(possible.fp),
                 possible.tp=sum(possible.tp),
                 errors=sum(fp+fn),
                 regions=length(tp),
                 row.names=NULL)
    })
  }#peaks.str
  error.totals <- do.call(rbind, error.by.peaks)
  exact <- converted$modelSelection
  exact$errors <- error.totals[paste(exact$peaks), "errors"]
  indices <- with(exact, {
    penaltyLearning::largestContinuousMinimumC(
      errors, max.log.lambda-min.log.lambda)
  })
  list(error.totals=error.totals,
       error.regions=error.regions.list,
       modelSelection=exact,
       target=c(exact$min.log.lambda[ indices[["start"]] ],
         exact$max.log.lambda[ indices[["end"]] ]))
### List of error.totals (data.frame with one row for each model size,
### with counts of incorrect labels), error.regions (list of
### data.frames with labels and error status for each model size),
### modelSelection (data.frame with one row for each model from
### exactModelSelection), target (numeric vector of length 2, lower
### and upper limits of target interval of log.lambda penalty values
### in the interval regression problem).
}, ex=function(){

  library(PeakSegJoint)
  data(H3K36me3.TDH.other.chunk1, envir=environment())
  lims <- c(43000000, 43200000) # left
  some.counts <-
    subset(H3K36me3.TDH.other.chunk1$counts,
           lims[1] < chromEnd & chromStart < lims[2])
  some.regions <-
    subset(H3K36me3.TDH.other.chunk1$regions,
           lims[1] < chromEnd & chromStart < lims[2])
  fit <- PeakSegJointSeveral(some.counts)
  converted <- ConvertModelList(fit)
  error.list <- PeakSegJointError(converted, some.regions)

  peaks.int.vec <- 1:3
  show.peaks <- subset(converted$peaks, peaks %in% peaks.int.vec)
  show.labels <- do.call(rbind, error.list$error.regions[paste(peaks.int.vec)])
  library(ggplot2)

  ann.colors <-
    c(noPeaks="#f6f4bf",
      peakStart="#ffafaf",
      peakEnd="#ff4c4c",
      peaks="#a445ee")
  ggplot()+
    penaltyLearning::geom_tallrect(aes(
      xmin=chromStart/1e3, xmax=chromEnd/1e3, fill=annotation),
                  alpha=0.5,
                  color="grey",
                  data=some.regions)+
    scale_fill_manual(values=ann.colors)+
    scale_linetype_manual("error type",
                          limits=c("correct", 
                                   "false negative",
                                   "false positive"
                                   ),
                          values=c(correct=0,
                                   "false negative"=3,
                                   "false positive"=1))+
    geom_step(aes(chromStart/1e3, count),
              color="grey50",
              data=some.counts)+
    penaltyLearning::geom_tallrect(aes(
      xmin=chromStart/1e3, xmax=chromEnd/1e3, linetype=status),
                  fill=NA,
                  color="black",
                  size=1,
                  data=show.labels)+
    geom_segment(aes(chromStart/1e3, 0,
                     xend=chromEnd/1e3, yend=0),
                 size=3,
                 color="deepskyblue",
                 data=show.peaks)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ peaks, scales="free")

})
