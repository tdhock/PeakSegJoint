PeakErrorSamples <- function
### Compute PeakError for several samples.
(peaks,
### data.frame of peaks with sample.id.
 regions
### data.frame of annotated region labels with sample.id.
 ){
  stopifnot(is.data.frame(peaks))
  stopifnot(is.data.frame(regions))
  getID <- function(df){
    with(df, paste0(sample.group, "/", sample.id))
  }
  regions.by.sample <- split(regions, getID(regions), drop=TRUE)
  peaks.by.sample <- split(peaks, getID(peaks), drop=TRUE)
  error.by.sample <- list()
  for(sample.path in names(regions.by.sample)){
    sample.regions <- regions.by.sample[[sample.path]]
    sample.peaks <- if(sample.path %in% names(peaks.by.sample)){
      peaks.by.sample[[sample.path]]
    }else{
      Peaks()
    }
    sample.error <- PeakErrorChrom(sample.peaks, sample.regions)
    if("chrom" %in% names(sample.regions)){
      sample.error$chrom <- sample.regions$chrom
    }
    sample.id <- sub(".*/", "", sample.path)
    sample.group <- sub("/.*", "", sample.path)
    error.by.sample[[sample.path]] <-
      data.frame(sample.id, sample.group, sample.error)
  }
  err <- do.call(rbind, error.by.sample)
  rownames(err) <- NULL
  err
### data.frame of error regions with sample.id.
}
