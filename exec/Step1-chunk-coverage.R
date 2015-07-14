library(data.table)

## Read one sample.bigwig file and save the coverage profile for
## each labeled chunk.

argv <-
  system.file(file.path("exampleData",
                        "bcell",
                        "McGill0322.bigwig"),
              package="PeakSegDP")

argv <- "~/exampleData/bcell/McGill0322.bigwig"

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

if(length(argv) != 1){
  stop("usage: Step1.R path/to/cellType/sample.bigwig")
}

bigwig.file <- argv[1]

bigwig.base <- basename(bigwig.file)
RData.base <- sub("bigwig$", "RData", bigwig.base)
type.dir <- dirname(bigwig.file)
cell.type <- basename(type.dir)
data.dir <- dirname(type.dir)
chunks.dir <- file.path(data.dir, "PeakSegJoint-chunks")
RData.path.vec <- Sys.glob(file.path(chunks.dir, "*", "regions.RData"))
for(RData.path in RData.path.vec){
  objs <- load(RData.path)
  stopifnot("chunk" %in% objs)
  tmp.bg <- tempfile()
  chunk.dir <- dirname(RData.path)
  out.RData <- file.path(chunk.dir, cell.type, RData.base)
  counts <-
    readBigWig(bigwig.file, chunk$chrom, chunk$chunkStart, chunk$chunkEnd)
  cat("Writing ", nrow(counts),
      " coverage rows to ", out.RData, "\n",
      sep="")
  out.dir <- dirname(out.RData)
  dir.create(out.dir, showWarnings=FALSE, recursive=TRUE)
  save(counts, file=out.RData)
}
cat("Done writing", length(RData.path.vec), "RData files.\n")
