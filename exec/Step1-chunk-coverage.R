if(!require(data.table))
  install.packages("data.table")

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
setkey(coverage, chrom, chromStart, chromEnd)

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
  counts <- foverlaps(coverage, chunk, nomatch=0L)
  setkey(counts, chromStart, chromEnd)
  chunk.dir <- dirname(RData.path)
  out.RData <- file.path(chunk.dir, cell.type, RData.base)
  cat("Writing ", nrow(counts),
      " coverage rows to ", out.RData, "\n",
      sep="")
  out.dir <- dirname(out.RData)
  dir.create(out.dir, showWarnings=FALSE, recursive=TRUE)
  save(counts, file=out.RData)
}
