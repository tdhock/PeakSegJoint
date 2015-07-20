library(PeakSegJoint)

argv <-
  c("~/exampleData/PeakSegJoint-chunks/trained.model.RData",
    "chr11")

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

PPN.cores()

if(length(argv) != 2){
  stop("usage: Step3.R PeakSegJoint-chunks/trained.model.RData chrom")
}

trained.model.RData <- argv[1]
chrom <- argv[2]

objs <- load(trained.model.RData)

chunks.dir <- dirname(trained.model.RData)
data.dir <- dirname(chunks.dir)
bigwig.file.vec <- Sys.glob(file.path(data.dir, "*", "*.bigwig"))

setkey(test.problems, chrom)
chrom.problems <- test.problems[chrom]

readBigWigSamples <- function(problem){
  counts.by.sample <- list()
  for(bigwig.file in bigwig.file.vec){
    sample.counts <- 
      readBigWig(bigwig.file, chrom,
                 problem$problemStart, problem$problemEnd)
    sample.id <- sub("[.]bigwig$", "", basename(bigwig.file))
    counts.by.sample[[sample.id]] <- with(sample.counts, {
      data.frame(chromStart, chromEnd, count)
    })
  }
  counts.by.sample
}

Step1Problem <- function(problem.i){
  ##cat(sprintf("%10d / %10d problems\n", problem.i, nrow(chrom.problems)))
  problem <- chrom.problems[problem.i, ]
  profile.list <- readBigWigSamples(problem)
  peak.only <- peak.or.null(profile.list)
  if(!is.null(peak.only)){
    data.table(problem, peak.only)
  }
}

step1.results.list <-
  mclapply.or.stop(seq_along(chrom.problems$problem.name), Step1Problem)

if(all(sapply(step1.results.list, is.null))){
  print(chrom.problems)
  warning("no computable models for any uniform size segmentation problems")
}else{
  step1.results <- do.call(rbind, step1.results.list)

  step2.problems <- clusterProblems(step1.results)

  ## ggplot()+
  ##   geom_segment(aes(problemStart, seq_along(problemStart),
  ##                    xend=problemEnd, yend=seq_along(problemStart)),
  ##                data=step1.results)+
  ##   geom_segment(aes(problemStart, -problem.i,
  ##                    xend=problemEnd, yend=-problem.i),
  ##                data=step2.problems)

  SegmentStep2 <- function(row.i){
    problem <- step2.problems[row.i, ]
    profile.list <- readBigWigSamples(problem)
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
      pred.log.lambda <- full.fit$predict(features)[[1]]
      selected <- 
        subset(info$modelSelection,
               min.log.lambda < pred.log.lambda &
               pred.log.lambda < max.log.lambda)
      stopifnot(nrow(selected) == 1)
      subset(info$peaks, peaks == selected$peaks)
    }
  }
  step2.peak.list <-
    mclapply.or.stop(seq_along(step2.problems$problem.name), SegmentStep2)

  pred.peaks <- do.call(rbind, step2.peak.list)

  if(0 == nrow(pred.peaks)){
    warning("no predicted peaks")
  }else{
    pred.peaks$chrom <- chrom
    
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

    sample.id.vec <- sort(unique(paste(pred.peaks$sample.id)))
    pred.peaks$peak.name <- with(pred.peaks, {
      sprintf("%s:%d-%d", chrom, chromStart, chromEnd)
    })
    peak.name.vec <- sort(unique(pred.peaks$peak.name))
    peak.mat <- 
      matrix(0, length(sample.id.vec), length(peak.name.vec),
             dimnames=list(sample.id=sample.id.vec,
               peak=peak.name.vec))
    i.mat <- with(pred.peaks, cbind(paste(sample.id), peak.name))
    peak.mat[i.mat] <- 1
    ## d.mat <- dist(peak.mat, method="manhattan")
    ## fit <- hclust(d.mat, method="average")
    ##plot(fit)

    ## combine chroms and save bed files later.
    ##peaks.by.sample <- split(pred.peaks, pred.peaks$sample.id)

    pred.dir <- file.path(data.dir, "PeakSegJoint-predictions")
    dir.create(pred.dir, showWarnings=FALSE)

    chrom.RData <- file.path(pred.dir, paste0(chrom, ".RData"))

    cat("Saved", ncol(peak.mat), "unique peaks across",
        nrow(peak.mat), "samples to", chrom.RData, "\n")

    save(pred.peaks, peak.mat, file=chrom.RData)
  }
}
