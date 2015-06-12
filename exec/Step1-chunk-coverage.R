library(data.table)

## Read one sample.bedGraph file and save the coverage profile for
## each labeled chunk.

argv <-
  system.file(file.path("exampleData",
                        "bcell",
                        "McGill0322.bedGraph"),
              package="PeakSegDP")

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

if(length(argv) != 1){
  stop("usage: Step1.R path/to/cellType/sample.bedGraph")
}

bedGraph.file <- argv[1]

coverage <- fread(bedGraph.file)
setnames(coverage, c("chrom", "chromStart", "chromEnd", "count"))
coverage[, chromStart1 := chromStart + 1L]
setkey(coverage, chrom, chromStart1, chromEnd)
chrom.ranges <-
  coverage[, .(min.chromStart=min(chromStart),
               max.chromEnd=max(chromEnd)),
           by=chrom]

bedGraph.base <- basename(bedGraph.file)
RData.base <- sub("bedGraph$", "RData", bedGraph.base)
type.dir <- dirname(bedGraph.file)
cell.type <- basename(type.dir)
data.dir <- dirname(type.dir)
chunks.dir <- file.path(data.dir, "PeakSegJoint-chunks")
RData.path.vec <- Sys.glob(file.path(chunks.dir, "*", "regions.RData"))
for(RData.path in RData.path.vec){
  objs <- load(RData.path)
  stopifnot("chunk" %in% objs)
  chunk[, chunkStart1 := chunkStart + 1L]
  setkey(chunk, chrom, chunkStart1, chunkEnd)
  counts.all <- foverlaps(coverage, chunk, nomatch=0L)
  counts <- counts.all[, .(chrom, chromStart, chromEnd, count)]
  chunk.dir <- dirname(RData.path)
  out.RData <- file.path(chunk.dir, cell.type, RData.base)
  cat("Writing ", nrow(counts),
      " coverage rows to ", out.RData, "\n",
      sep="")
  out.dir <- dirname(out.RData)
  dir.create(out.dir, showWarnings=FALSE, recursive=TRUE)
  save(counts, chrom.ranges, file=out.RData)
}
