library(PeakSegJoint)

argv <- "~/genomepipelines/H3K4me3_TDH_immune/PeakSegJoint-predictions"

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

if(length(argv) != 1){
  stop("usage: Step6.R path/to/PeakSegJoint-predictions")
}

pred.dir <- normalizePath(argv, mustWork=TRUE)

RData.file.vec <- Sys.glob(file.path(pred.dir, "*.RData"))

peak.mat.list <- list()
pred.peak.list <- list()
for(RData.file in RData.file.vec){
  objs <- load(RData.file)
  pred.peak.list[[RData.file]] <- data.frame(pred.peaks)
  peak.mat.list[[RData.file]] <- peak.mat
}
all.peaks.mat <- do.call(cbind, peak.mat.list)
all.peaks.df <- do.call(rbind, pred.peak.list)
all.peaks.df$sample.path <-
  with(all.peaks.df, file.path(sample.group, sample.id))

data.dir <- dirname(pred.dir)

predictions.RData <- file.path(data.dir, "PeakSegJoint.predictions.RData")
predictions.csv <- file.path(data.dir, "PeakSegJoint.predictions.csv")
save(all.peaks.mat, all.peaks.df, file=predictions.RData)
write.csv(all.peaks.mat, predictions.csv, quote=TRUE,
          row.names=TRUE)

bigwig.file.vec <- Sys.glob(file.path(data.dir, "*", "*.bigwig"))
chrom.ranges <- bigWigInfo(bigwig.file.vec[1])
just.ends <- data.frame(chrom.ranges)[, c("chrom", "chromEnd")]
chrom.file <- tempfile()
write.table(just.ends, chrom.file, quote=FALSE,
            row.names=FALSE, col.names=FALSE)

peaks.by.path <- split(all.peaks.df, all.peaks.df$sample.path)
for(sample.path in names(peaks.by.sample)){
  sample.peaks <- peaks.by.path[[sample.path]]
  ord <- with(sample.peaks, order(chrom, chromStart))
  sorted.peaks <- sample.peaks[ord, ]
  bed.peaks <- sorted.peaks[, c("chrom", "chromStart", "chromEnd")]
  rownames(bed.peaks) <- NULL
  bed.file <- file.path(data.dir, paste0(sample.path, ".bed"))
  bigBed.file <- file.path(data.dir, paste0(sample.path, ".bigBed"))
  write.table(bed.peaks, bed.file, quote=FALSE, row.names=FALSE, col.names=FALSE)
  s <- ifelse(nrow(bed.peaks)==1, "", "s")
  message("wrote ", nrow(bed.peaks), " peak", s, " to", bed.file)
  bigBedCmd <- paste("bedToBigBed", bed.file, chrom.file, bigBed.file)
  system(bigBedCmd)
}
