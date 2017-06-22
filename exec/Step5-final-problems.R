library(PeakSegJoint)
library(data.table)

argv <-
  c("~/exampleData/combined.problems.RData",
    "1")

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

ppn <- PPN.cores()
if(!is.na(ppn))options(mc.cores=ppn/2)

if(length(argv) != 2){
  stop("usage: Step5.R PeakSegJoint-chunks/combined.problems.RData jobNum")
}

combined.problems.RData <- normalizePath(argv[1], mustWork=TRUE)
jobNum <- argv[2]

objs <- load(combined.problems.RData)
job.problems <- combined.problems[[jobNum]]

data.dir <- dirname(combined.problems.RData)
chunks.dir <- file.path(data.dir, "PeakSegJoint-chunks")
trained.model.RData <- file.path(chunks.dir, "trained.model.RData")
load(trained.model.RData)
bases.per.problem <- train.errors.picked$bases.per.problem

bigwig.file.vec <- Sys.glob(file.path(data.dir, "*", "*.bigwig"))
sample.group.vec <- basename(dirname(bigwig.file.vec))
sample.id.vec <- sub("[.]bigwig$", "", basename(bigwig.file.vec))
names(bigwig.file.vec) <- paste0(sample.group.vec, "/", sample.id.vec)

SegmentFinal <- function(row.i){
  message(sprintf("%10d / %10d final problems",
                  row.i, nrow(job.problems)))
  problem <- job.problems[row.i, ]
  profile.list <- readBigWigSamples(problem, bigwig.file.vec)
  fit <- tryCatch({
    PeakSegJointSeveral(profile.list)
  }, error=function(e){
    NULL
  })
  if(is.null(fit)){
    return(NULL)
  }
  info <- ConvertModelList(fit)
  if(is.data.frame(info$peaks)){
    features <- featureMatrix(profile.list)
    fmat <- rbind(colSums(features))
    pred.log.lambda <- full.fit$predict(fmat)
    selected <- 
      subset(info$modelSelection,
             min.log.lambda < pred.log.lambda &
               pred.log.lambda < max.log.lambda)
    stopifnot(nrow(selected) == 1)
    peak.df <- subset(info$peaks, peaks == selected$peaks)
    if(nrow(peak.df)){
      data.table(problem, peak.df)
    }
  }
}
final.peak.list <-
  mclapply.or.stop(seq_along(job.problems$problem.name), SegmentFinal)
pred.peaks <- do.call(rbind, final.peak.list)

if(is.null(pred.peaks) || 0 == nrow(pred.peaks)){
  warning("no predicted peaks")
}else{
  ## If for some reason we do not have sample.group info, there is a
  ## problem!
  stopifnot(all(!is.na(pred.peaks$sample.group)))

  ## big.problem <- with(step1.results, {
  ##   data.table(problemStart=min(problemStart),
  ##              problemEnd=max(problemEnd))
  ## })
  ## all.counts.list <- readBigWigSamples(big.problem)
  ## for(sample.id in names(all.counts.list)){
  ##   all.counts.list[[sample.id]]$sample.id <- sample.id
  ## }
  ## all.counts <- data.table(do.call(rbind, all.counts.list))
  ## ggplot()+
  ##   geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
  ##                 ymin=0, ymax=count),
  ##             data=all.counts,
  ##             color="grey")+
  ##   geom_segment(aes(chromStart/1e3, 0,
  ##                    xend=chromEnd/1e3, yend=0),
  ##                data=pred.peaks,
  ##                color="deepskyblue",
  ##                size=2)+
  ##   theme_bw()+
  ##   theme(panel.margin=grid::unit(0, "cm"))+
  ##   facet_grid(sample.id ~ ., scales="free")

  sample.id.vec <- names(bigwig.file.vec)
  pred.peaks$peak.name <- with(pred.peaks, {
    sprintf("%s:%d-%d", chrom, chromStart, chromEnd)
  })
  peak.name.vec <- sort(unique(pred.peaks$peak.name))
  peak.mat <- 
    matrix(0, length(sample.id.vec), length(peak.name.vec),
           dimnames=list(sample.id=sample.id.vec,
             peak=peak.name.vec))
  i.mat <- with(pred.peaks, {
    cbind(paste0(sample.group, "/", sample.id), peak.name)
  })
  peak.mat[i.mat] <- 1
  ## d.mat <- dist(peak.mat, method="manhattan")
  ## fit <- hclust(d.mat, method="average")
  ##plot(fit)

  ## combine chroms and save bed files later.
  ##peaks.by.sample <- split(pred.peaks, pred.peaks$sample.id)

  pred.dir <- file.path(data.dir, "PeakSegJoint-predictions")

  job.RData <- file.path(pred.dir, paste0(jobNum, ".RData"))

  cat("Saved", ncol(peak.mat), "unique peaks across",
      nrow(peak.mat), "samples to", job.RData, "\n")

  save(pred.peaks, peak.mat, file=job.RData)
}
