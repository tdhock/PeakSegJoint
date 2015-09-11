PeakSegJointError <- function
### Compute number of incorrect regions for every PeakSegJoint model.
(converted,
### Result of ConvertModelList.
 problem.regions
### data.frame of annotated region labels.
 ){
  getID <- function(df){
    with(df, paste0(sample.group, "/", sample.id))
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
        data.frame(sample.id, sample.group, model.error)
    }
    model.error <- do.call(rbind, error.by.sample)
    error.regions.list[[peaks.str]] <- model.error
    error.by.peaks[[peaks.str]] <- with(model.error, {
      data.frame(peaks=as.integer(peaks.str),
                 fp=sum(fp),
                 fn=sum(fn),
                 tp=sum(tp),
                 possible.fp=sum(possible.fp),
                 possible.tp=sum(possible.tp),
                 errors=sum(fp+fn),
                 regions=length(tp))
    })
  }#peaks.str
  error.totals <- do.call(rbind, error.by.peaks)
  exact <- converted$modelSelection
  exact$errors <- error.totals[paste(exact$peaks), "errors"]
  indices <- with(exact, {
    largestContinuousMinimum(errors, max.log.lambda-min.log.lambda)
  })
  list(error.totals=error.totals,
       error.regions=error.regions.list,
       modelSelection=exact,
       target=c(exact$min.log.lambda[indices$start],
         exact$max.log.lambda[indices$end]))
### List of error.totals, error.regions, modelSelection, target.
}
