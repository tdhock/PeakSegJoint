library(PeakSegJoint)

argv <- system.file("exampleData", "PeakSegJoint-predictions",
                    mustWork=TRUE,
                    package="PeakSegJoint")

argv <- "~/genomepipelines/H3K4me3_TDH_immune/PeakSegJoint-predictions"

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

if(length(argv) != 1){
  stop("usage: Step6.R path/to/PeakSegJoint-predictions")
}

pred.dir <- normalizePath(argv[1], mustWork=TRUE)
data.dir <- dirname(pred.dir)

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
rownames(all.peaks.df) <- NULL
all.peaks.df$sample.path <-
  with(all.peaks.df, file.path(sample.group, sample.id))

writeBoth <- function(peaks.df, bed.path){
  bigBed.file <- sub("[.]bed$", ".bigBed", bed.path)
  write.table(peaks.df, bed.path, quote=FALSE, row.names=FALSE, col.names=FALSE)
  s <- ifelse(nrow(peaks.df)==1, "", "s")
  message("wrote ", nrow(peaks.df), " peak", s, " to ", bed.path)
  bigBedCmd <- paste("bedToBigBed", bed.path, chrom.file, bigBed.file)
  system(bigBedCmd)
}
bigwig.file.vec <- Sys.glob(file.path(data.dir, "*", "*.bigwig"))
chrom.ranges <- bigWigInfo(bigwig.file.vec[1])
just.ends <- data.frame(chrom.ranges)[, c("chrom", "chromEnd")]
chrom.file <- tempfile()
write.table(just.ends, chrom.file, quote=FALSE,
            row.names=FALSE, col.names=FALSE)

u.cols <- c("chrom", "chromStart", "chromEnd")
u.peaks <-
  unique(all.peaks.df[, u.cols])
rownames(u.peaks) <- with(u.peaks, {
  paste0(chrom, ":", chromStart, "-", chromEnd)
})
count.mat <- with(all.peaks.df, table(sample.group, peak.name))
getCount <- function(x){
  count.str <- paste0(x, ":", rownames(count.mat))
  non.zero <- count.str[x != 0]
  paste(non.zero, collapse="/")
}
str.vec <- apply(count.mat, 2, getCount)
u.peaks$name <- count.vec[rownames(u.peaks)]
count.tab <- colSums(all.peaks.mat)
u.peaks$score <- count.tab[rownames(u.peaks)]
summary.bed <- file.path(data.dir, "PeakSegJoint.summary.bed")
writeBoth(u.peaks, summary.bed)

predictions.RData <- file.path(data.dir, "PeakSegJoint.predictions.RData")
predictions.csv <- file.path(data.dir, "PeakSegJoint.predictions.csv")
save(all.peaks.mat, all.peaks.df, file=predictions.RData)
write.csv(all.peaks.mat, predictions.csv, quote=TRUE,
          row.names=TRUE)

peaks.by.path <- split(all.peaks.df, all.peaks.df$sample.path)
for(sample.path in names(peaks.by.path)){
  sample.peaks <- peaks.by.path[[sample.path]]
  ord <- with(sample.peaks, order(chrom, chromStart))
  sorted.peaks <- sample.peaks[ord, ]
  bed.peaks <- sorted.peaks[, c("chrom", "chromStart", "chromEnd")]
  rownames(bed.peaks) <- NULL
  bed.file <- file.path(data.dir, paste0(sample.path, ".bed"))
  writeBoth(bed.peaks, bed.file)
}
