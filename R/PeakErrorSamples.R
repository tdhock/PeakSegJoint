PeakErrorSamples <- function
### Compute PeakError for several samples.
(peaks,
### data.frame of peaks with sample.id.
 regions
### data.frame of annotated region labels with sample.id.
 ){
  stopifnot(is.data.frame(peaks))
  stopifnot(is.data.frame(regions))
  regions.by.sample <- split(regions, regions$sample.id, drop=TRUE)
  peaks.by.sample <- split(peaks, peaks$sample.id, drop=TRUE)
  error.by.sample <- list()
  for(sample.id in names(regions.by.sample)){
    sample.regions <- regions.by.sample[[sample.id]]
    sample.peaks <- if(sample.id %in% names(peaks.by.sample)){
      peaks.by.sample[[sample.id]]
    }else{
      Peaks()
    }
    sample.error <- PeakErrorChrom(sample.peaks, sample.regions)
    error.by.sample[[sample.id]] <-
      data.frame(sample.id, sample.error)
  }
  err <- do.call(rbind, error.by.sample)
  rownames(err) <- NULL
  err
### data.frame of error regions with sample.id.
}
