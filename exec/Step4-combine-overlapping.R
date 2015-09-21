library(PeakSegJoint)

argv <-
  c("~/exampleData/PeakSegJoint-overlapping", "200")

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

PPN.cores()

if(length(argv) != 2){
  stop("usage: Step4.R path/to/PeakSegJoint-overlapping numJobs")
}

overlapping.dir <- normalizePath(argv[1], mustWork=TRUE)
numJobs <- as.integer(argv[2])

data.dir <- dirname(overlapping.dir)

job.RData.vec <- Sys.glob(file.path(overlapping.dir, "*.RData"))
peaks.by.job <- list()
for(job.RData in job.RData.vec){
  objs <- load(job.RData)
  peaks.by.job[[job.RData]] <- overlapping.peaks
}
all.overlapping.peaks <- do.call(rbind, peaks.by.job)
message(nrow(all.overlapping.peaks), " peaks in ",
        length(job.RData.vec), " RData files.")

peaks.by.chrom <- split(all.overlapping.peaks, all.overlapping.peaks$chrom)
clustered.by.chrom <- list()
for(chrom in names(peaks.by.chrom)){
  chrom.peaks <- peaks.by.chrom[[chrom]]
  clustered.peaks <- clusterProblems(chrom.peaks)
  clustered.by.chrom[[chrom]] <-
    data.table(chrom, clustered.peaks)
}
final.problems <- do.call(rbind, clustered.by.chrom)
final.problems$job <- sort(rep(1:numJobs, l=nrow(final.problems)))

message(nrow(final.problems), " final segmentation problems.")

combined.problems <- split(final.problems, final.problems$job)

out.RData <- file.path(data.dir, "combined.problems.RData")

save(combined.problems, file=out.RData)
