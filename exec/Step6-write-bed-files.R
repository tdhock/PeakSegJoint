library(PeakSegJoint)

argv <- "~/genomepipelines/H3K4me3_TDH_immune/PeakSegJoint-predictions"

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

pred.dir <- normalizePath(argv, mustWork=TRUE)

RData.file.vec <- Sys.glob(file.path(pred.dir, "*.RData"))

peak.mat.list <- list()
pred.peak.list <- list()
for(RData.file in RData.file.vec){
  objs <- load(RData.file)
  pred.peak.list[[RData.file]] <- pred.peaks
  peak.mat.list[[RData.file]] <- peak.mat
}
all.peaks.mat <- do.call(cbind, peak.mat.list)
all.peaks.df <- do.call(rbind, pred.peak.list)

data.dir <- dirname(pred.dir)

predictions.RData <- file.path(data.dir, "PeakSegJoint.predictions.RData")
predictions.csv <- file.path(data.dir, "PeakSegJoint.predictions.csv")
save(all.peaks.mat, all.peaks.df, file=predictions.RData)
write.csv(all.peaks.mat, predictions.csv, quote=TRUE,
          row.names=TRUE)

bigwig.file.vec <- Sys.glob(file.path(data.dir, "*", "*.bigwig"))
chrom.ranges <- bigWigInfo(bigwig.file.vec[1])
chrom.file <- tempfile()
write.table(chrom.ranges[, c("chrom", "chromEnd")], chrom.file, quote=FALSE,
            row.names=FALSE, col.names=FALSE)

peaks.by.sample <- split(all.peaks.df, all.peaks.df$sample.id)
for(sample.id in names(peaks.by.sample)){
  sample.peaks <- peaks.by.sample[[sample.id]]
  bed.peaks <- sample.peaks[, c("chrom", "chromStart", "chromEnd")]
  rownames(bed.peaks) <- NULL
  bed.file <- file.path(pred.dir, paste0(sample.id, ".bed"))
  bigBed.file <- file.path(pred.dir, paste0(sample.id, ".bigBed"))
  write.table(bed.peaks, bed.file, quote=FALSE, row.names=FALSE, col.names=FALSE)
  bigBedCmd <- paste("bedToBigBed", bed.file, chrom.file, bigBed.file)
}
