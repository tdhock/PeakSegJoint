library(PeakSegJoint)
library(PeakError)

## Compute PeakSegJoint segmentations for one labeled chunk in the
## train data set.

argv <-
  system.file(file.path("exampleData",
                        "PeakSegJoint-chunks",
                        "1"),
              package="PeakSegJoint")

argv <- "~/exampleData/PeakSegJoint-chunks/3"
argv <- "~/projects/PeakSegJoint-paper/PeakSegJoint-chunks/H3K36me3_AM_immune/21"
argv <- "~/projects/PeakSegJoint-paper/PeakSegJoint-chunks/H3K4me3_PGP_immune/2"

argv <- commandArgs(trailingOnly=TRUE)

ppn <- PPN.cores()
if(!is.na(ppn))options(mc.cores=ppn/2)

print(argv)

if(length(argv) != 1){
  stop("usage: Step1.R path/to/PeakSegJoint-chunks/012354")
}

chunk.dir <- argv[1]
bins.per.problem <- 500L
chunk.id <- basename(chunk.dir)
chunks.dir <- dirname(chunk.dir)
data.dir <- dirname(chunks.dir)
bigwig.file.vec <- Sys.glob(file.path(data.dir, "*", "*.bigwig"))

regions.RData <- file.path(chunk.dir, "regions.RData")

objs <- load(regions.RData)
regions$region.i <- 1:nrow(regions)
chrom <- paste(regions$chrom[1])
regions[, chromStart1 := chromStart + 1L]
regions[, id.group := paste(sample.id, sample.group)]

counts.by.sample <- list()
for(bigwig.file in bigwig.file.vec){
  sample.counts <-
    readBigWig(bigwig.file, chunk$chrom, chunk$chunkStart, chunk$chunkEnd)
  sample.id <- sub("[.]bigwig$", "", basename(bigwig.file))
  sample.group <- basename(dirname(bigwig.file))
  counts.by.sample[[paste(sample.id, sample.group)]] <-
    data.table(sample.id, sample.group, sample.counts)
}
counts <- do.call(rbind, counts.by.sample)
counts[, chromStart1 := chromStart + 1L]
setkey(counts, chromStart1, chromEnd)

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")
sample.problems.dt <- do.call(rbind, problems.by.res)

## ggplot()+
##   ggtitle(regions.RData)+
##   geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3, fill=annotation),
##                 color="grey",
##                 alpha=0.5,
##                 data=regions)+
##   geom_segment(aes(problemStart/1e3, seq_along(bases.per.problem),
##                    xend=problemEnd/1e3, yend=seq_along(bases.per.problem)),
##                data=data.table(sample.problems.dt,
##                                sample.group="problems", sample.id="problems"))+
##   scale_fill_manual(values=ann.colors)+
##   geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
##                 ymin=0, ymax=count),
##             fill="grey50",
##             data=counts)+
##   theme_bw()+
##   theme(panel.margin=grid::unit(0, "cm"))+
##   facet_grid(sample.group + sample.id ~ ., scales="free")

SampleProblems <- function(problem.i){
  cat(sprintf("%4d / %4d problems\n", problem.i, nrow(sample.problems.dt)))
  problem <- sample.problems.dt[problem.i, ]
  bases.per.bin <- as.integer(problem$bases.per.problem/bins.per.problem)
  problem.name <- paste(problem$problem.name)
  problem.regions <-
    regions[! (chromEnd < problem$problemStart |
                 problem$problemEnd < chromStart), ]
  setkey(problem.regions, id.group)
  
  ## ggplot()+
  ##   ggtitle(problem.name)+
  ##   geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3, fill=annotation),
  ##                 color="grey",
  ##                 alpha=0.5,
  ##                 data=regions)+
  ##   geom_segment(aes(problemStart/1e3, 0,
  ##                    xend=problemEnd/1e3, yend=0),
  ##                data=problem)+
  ##   scale_fill_manual(values=ann.colors)+
  ##   geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
  ##                 ymin=0, ymax=count),
  ##             fill="grey50",
  ##             data=counts)+
  ##   theme_bw()+
  ##   theme(panel.margin=grid::unit(0, "cm"))+
  ##   facet_grid(sample.group + sample.id ~ ., scales="free")

  models.by.sample <- list()
  for(id.group in names(counts.by.sample)){
    sample.counts <- counts.by.sample[[id.group]]
    problem.counts <-
      sample.counts[! (chromEnd < problem$problemStart |
                         problem$problemEnd < chromStart), ]
    start <- as.integer(problem$problemStart)
    model <- segmentBins(
      problem.counts, start, bases.per.bin, bins.per.problem)
    sample.regions <- problem.regions[id.group]
    if(!is.na(sample.regions$sample.id[1])){
      
      ## ggplot()+
      ##   ggtitle(problem.name)+
      ##   geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
      ##                     fill=annotation),
      ##                 color="grey",
      ##                 alpha=0.5,
      ##                 data=sample.regions)+
      ##   geom_segment(aes(problemStart/1e3, -1,
      ##                    xend=problemEnd/1e3, yend=-1),
      ##                size=2,
      ##                data=problem)+
      ##   scale_fill_manual(values=ann.colors)+
      ##   geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
      ##                 ymin=0, ymax=count),
      ##             fill="grey50",
      ##             data=problem.counts)+
      ##   geom_line(aes((chromStart+chromEnd)/2e3, mean),
      ##             data=model$bins)+
      ##   theme_bw()+
      ##   theme(panel.margin=grid::unit(0, "cm"))+
      ##   facet_grid(sample.group + sample.id ~ ., scales="free")
      
      model$fit$error$incorrect.regions <- NA
      for(peaks.str in names(model$fit$peaks)){
        peaks <- model$fit$peaks[[peaks.str]]
        error.regions <- PeakErrorChrom(peaks, sample.regions)
        model$fit$error[peaks.str, "incorrect.regions"] <-
          with(error.regions, sum(fp+fn))
      }
      model$modelSelection$incorrect.regions <-
        model$fit$error[paste(model$modelSelection$peaks), "incorrect.regions"]
      target.indices <- with(model$modelSelection, largestContinuousMinimum(
        incorrect.regions, max.log.lambda-min.log.lambda))
      model$target <- with(model$modelSelection, c(
        min.log.lambda[target.indices$start],
        max.log.lambda[target.indices$end]))
    }
    models.by.sample[[id.group]] <- model
  }#id.group
  models.by.sample
}#problem.i
models.by.problem <-
  mclapply.or.stop(seq_along(sample.problems.dt$problem.name), SampleProblems)

out.RData <- file.path(chunk.dir, "models.by.problem.RData")
save(models.by.problem, sample.problems.dt, file=out.RData)
