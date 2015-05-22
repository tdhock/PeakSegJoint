PeakSegJointError <- function
### Compute number of incorrect regions for every PeakSegJoint model.
(converted,
### Result of ConvertModelList.
 problem.regions
### data.frame of annotated region labels.
 ){
  regions.by.sample <- split(problem.regions, problem.regions$sample.id)
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
      split(model.peaks, model.peaks$sample.id)
    }
    error.by.sample <- list()
    for(sample.id in names(regions.by.sample)){
      sample.regions <- regions.by.sample[[sample.id]]
      sample.peaks <- if(sample.id %in% names(peaks.by.sample)){
        peaks.by.sample[[sample.id]]
      }else{
        Peaks()
      }
      model.error <- PeakErrorChrom(sample.peaks, sample.regions)
      error.by.sample[[sample.id]] <-
        data.frame(sample.id, model.error)
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
