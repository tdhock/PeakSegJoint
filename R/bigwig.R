readBigWig <- function
### Read part of a bigWig file into R as a data.table (assumes
### bigWigToBedGraph is present on your PATH).
(bigwig.file,
### path or URL of bigwig file.
 chrom,
### chromosome to read.
 start,
### position before reading.
 end,
### last position to read.
 bedGraph.file=tempfile()
### plain text file where coverage is saved before reading into R.
 ){
  stopifnot(length(bigwig.file) == 1)
  stopifnot(length(chrom) == 1)
  stopifnot(length(start) == 1)
  stopifnot(length(end) == 1)
  stopifnot(length(bedGraph.file) == 1)
  stopifnot(is.character(bigwig.file))
  stopifnot(is.character(bedGraph.file))
  stopifnot(is.character(chrom))
  start <- as.integer(start)
  end <- as.integer(end)
  stopifnot(0 <= start)
  stopifnot(start < end)
  stopifnot(end < Inf)
  cmd <-
    sprintf("bigWigToBedGraph -chrom=%s -start=%d -end=%d %s %s",
            chrom, start, end,
            bigwig.file, bedGraph.file)
  status <- system(cmd)
  if(status != 0){
    stop("error code ", status, " for\n", cmd)
  }
  if(file.info(bedGraph.file)$size == 0){
    data.table(chrom=character(),
               chromStart=integer(),
               chromEnd=integer(),
               count=integer())
  }else{
    bg <- fread(bedGraph.file)
    setnames(bg, c("chrom", "chromStart", "chromEnd", "norm"))
    stopifnot(0 <= bg$norm)
    nonzero <- bg[0 < norm, ]
    min.nonzero.norm <- min(nonzero[, norm])
    nonzero[, count := as.integer(norm/min.nonzero.norm) ]
    nonzero[, .(
      chrom,
      chromStart,
      chromEnd,
      count
      )]
  }
### data.table with columns chrom chromStart chromEnd count.
}

bigWigInfo <- function
### Run bigWigInfo to find chrom sizes.
(bigwig.file
### path or URL of bigwig file.
 ){
  stopifnot(is.character(bigwig.file))
  stopifnot(length(bigwig.file) == 1)
  cmd <- paste("bigWigInfo", bigwig.file, "-chroms")
  info.lines <- system(cmd, intern=TRUE)
  chrom.lines <- grep("^\t", info.lines, value=TRUE)
  chrom.df <- read.table(text=chrom.lines)
  names(chrom.df) <- c("chrom", "chromStart", "chromEnd")
  chrom.df$chrom <- paste(chrom.df$chrom)
  chrom.df
}

