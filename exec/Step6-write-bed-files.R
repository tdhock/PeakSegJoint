library(PeakSegJoint)
library(data.table)

argv <- system.file("exampleData", "PeakSegJoint-predictions",
                    package="PeakSegJoint")

argv <- "~/genomepipelines/H3K4me3_TDH_immune/PeakSegJoint-predictions"

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

if(length(argv) != 1){
  stop("usage: Step6.R path/to/PeakSegJoint-predictions")
}

pred.dir <- normalizePath(argv[1], mustWork=TRUE)
data.dir <- dirname(pred.dir)

positive.regions.RData <- file.path(data.dir, "positive.regions.RData")
load(positive.regions.RData)

RData.file.vec <- Sys.glob(file.path(pred.dir, "*.RData"))

peak.mat.list <- list()
pred.peak.list <- list()
for(RData.file in RData.file.vec){
  objs <- load(RData.file)
  pred.peak.list[[RData.file]] <- data.table(pred.peaks)
  peak.mat.list[[RData.file]] <- peak.mat
}
all.peaks.mat <- do.call(cbind, peak.mat.list)
all.peaks.dt <- do.call(rbind, pred.peak.list)
all.peaks.dt[, sample.path := file.path(sample.group, sample.id)]
all.peaks.df <- data.frame(all.peaks.dt)

## Join positive region labels with peak calls.
Input.counts <- all.peaks.dt[, list(
  Input=sum(sample.group=="Input")
  ), by=.(peak.name, chrom, chromStart, chromEnd)]
if(is.data.frame(positive.regions)){
  setkey(Input.counts, chrom, chromStart, chromEnd)
  positive.dt <- data.table(positive.regions)
  setkey(positive.dt, chrom, regionStart, regionEnd)
  over.dt <- foverlaps(Input.counts, positive.dt, nomatch=0L)
  thresh.dt <- over.dt[, list(
    specific=sum(!input.has.peak),
    nonspecific=sum(input.has.peak)
    ), by=Input]
  ## The decision function calls a peak "specific" if at most T input
  ## samples have peaks, where T is the learned threshold: in the error
  ## rate table below, [min.thresh, max.thresh), e.g. [4, 0) means that
  ## thresholds 4, 3, 2, 1 all have the same error. 
  specific.error.dt <- thresh.dt[, {
    min.thresh <- c(-Inf, Input)
    max.thresh <- c(Input, Inf)
    threshold <- floor((min.thresh+max.thresh)/2)
    threshold[1] <- min(Input)
    threshold[length(threshold)] <- max(Input)
    data.table(
      min.thresh,
      max.thresh,
      threshold,
      FN=c(0, cumsum(rev(specific))),
      FP=c(rev(cumsum(nonspecific)), 0))
  }]
  specific.error.dt[, errors := FP + FN]
  max.input.samples <- specific.error.dt[which.min(errors), threshold]
  specific.error <- data.frame(specific.error.dt)
  print(specific.error)
  cat('A peak is called "specific" and included in output files if at most',
      max.input.samples,
      "Input samples have peaks. Counts:\n")
  Input.counts[, call := ifelse(Input <= max.input.samples,
                        "specific", "nonspecific")]
  print(table(Input.counts$call))
  specific.peak.names <- Input.counts[call=="specific", peak.name]
}else{
  cat("All peaks included in output files,",
      "since there are no labeled peaks in Input samples.\n")
  specific.peak.names <- colnames(all.peaks.mat)
  specific.error <- NULL
}

writeBoth <- function(peaks.df, bed.path, header.line){
  bigBed.file <- sub("[.]bed$", ".bigBed", bed.path)
  write.table(peaks.df, bed.path, 
              quote=FALSE, row.names=FALSE, col.names=FALSE)
  s <- ifelse(nrow(peaks.df)==1, "", "s")
  message("wrote ", nrow(peaks.df), " peak", s, " to ", bed.path)
  bigBedCmd <- paste("bedToBigBed", bed.path, chrom.file, bigBed.file)
  status <- system(bigBedCmd)
  if(status != 0){
    stop("error in bedToBigBed")
  }
  ## Also write a bed.gz file with a header for uploading as a custom
  ## track.
  bed.gz <- paste0(bed.path, ".gz")
  con <- gzfile(bed.gz, "w")
  writeLines(header.line, con)
  write.table(peaks.df, con, 
              quote=FALSE, row.names=FALSE, col.names=FALSE)
  close(con)
}
bigwig.file.vec <- Sys.glob(file.path(data.dir, "*", "*.bigwig"))
chrom.ranges <- bigWigInfo(bigwig.file.vec[1])
just.ends <- data.frame(chrom.ranges)[, c("chrom", "chromEnd")]
chrom.file <- tempfile()
write.table(just.ends, chrom.file, quote=FALSE,
            row.names=FALSE, col.names=FALSE)

## Save a summary bed file for easy viewing of the entire model as a
## custom track on UCSC.
specific.peaks.df <- subset(all.peaks.df, peak.name %in% specific.peak.names)
u.cols <- c("chrom", "chromStart", "chromEnd")
u.peaks <-
  unique(specific.peaks.df[, u.cols])
rownames(u.peaks) <- with(u.peaks, {
  paste0(chrom, ":", chromStart, "-", chromEnd)
})
specific.peaks.df$group <- sub("[0-9]*$", "", specific.peaks.df$sample.group)
count.mat <- with(specific.peaks.df, table(group, peak.name))
getCount <- function(x){
  count.str <- paste0(x, ":", rownames(count.mat))
  non.zero <- count.str[x != 0]
  paste(non.zero, collapse="/")
}
str.vec <- apply(count.mat, 2, getCount)
u.peaks$name <- str.vec[rownames(u.peaks)]
count.tab <- colSums(all.peaks.mat)
u.peaks$score <- count.tab[rownames(u.peaks)]
summary.bed <- file.path(data.dir, "PeakSegJoint.summary.bed")
header <- 
  paste("track",
        "visibility=pack",
        "name=PeakSegJointGroups",
        'description="PeakSegJoint sample counts by group"')
writeBoth(u.peaks, summary.bed, header)

predictions.RData <- file.path(data.dir, "PeakSegJoint.predictions.RData")
predictions.csv <- file.path(data.dir, "PeakSegJoint.predictions.csv")
peak.Input.counts <- data.frame(Input.counts)
save(all.peaks.mat, all.peaks.df,
     specific.peak.names,
     peak.Input.counts, specific.error,
     file=predictions.RData)
specific.peaks.mat <- all.peaks.mat[, specific.peak.names]
write.csv(specific.peaks.mat, predictions.csv, quote=TRUE,
          row.names=TRUE)

peaks.by.path <- split(specific.peaks.df, specific.peaks.df$sample.path)
for(sample.path in names(peaks.by.path)){
  sample.peaks <- peaks.by.path[[sample.path]]
  ord <- with(sample.peaks, order(chrom, chromStart))
  sorted.peaks <- sample.peaks[ord, ]
  bed.peaks <- sorted.peaks[, c("chrom", "chromStart", "chromEnd")]
  rownames(bed.peaks) <- NULL
  bed.file <- file.path(data.dir, paste0(sample.path, ".bed"))
  desc.line <-
    sprintf('description="PeakSegJoint calls for sample %s"',
            sample.path)
  header <- 
    paste("track",
          "visibility=pack",
          paste0("name=PeakSegJoint:", sample.path),
          desc.line)
  writeBoth(bed.peaks, bed.file, header)
}
