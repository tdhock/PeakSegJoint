library(ggplot2)
library(xtable)
library(PeakSegJoint)
options(xtable.print.results=FALSE)

argv <-
  c(system.file(file.path("exampleData",
                          "PeakSegJoint-chunks"),
                package="PeakSegJoint"),
    "200")

argv <- c("PeakSegJoint-chunks/H3K36me3_TDH_immune", "200")
argv <- c("~/exampleData/PeakSegJoint-chunks", "200")
argv <- c("~/projects/H3K27ac_TDH/PeakSegJoint-chunks", "200")

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

PPN.cores()

if(length(argv) != 2){
  stop("usage: Step2.R path/to/PeakSegJoint-chunks numJobs")
}

chunks.dir <- normalizePath(argv[1])
numJobs <- as.integer(argv[2])
data.dir <- dirname(chunks.dir)

problems.RData.vec <- Sys.glob(file.path(chunks.dir, "*", "problems.RData"))

err.mat.list <- list()
chunk.best.list <- list()
bpp.list <- list()
for(chunk.i in seq_along(problems.RData.vec)){
  problems.RData <- problems.RData.vec[[chunk.i]]
  chunk.dir <- dirname(problems.RData)
  chunk.id <- basename(chunk.dir)
  index.html <-
    file.path("..", chunk.id, "figure-train-errors", "index.html")
  first.png <-
    file.path("..", chunk.id, "figure-train-errors.png")
  objs <- load(problems.RData)
  err.vec <- res.error$errors
  min.df <- subset(res.error, errors==min(errors))
  best.errors <- min(res.error$errors)
  bpp.list[[chunk.id]] <- min.df$bases.per.problem
  href <- sprintf('<a href="%s">
  <img src="%s" alt="chunk%s" />
</a>', index.html, first.png, chunk.id)
  chunk.best.list[[chunk.id]] <-
    data.frame(selected.model=href,
               best.errors,
               regions=res.error$regions[1])
  names(err.vec) <- res.error$bases.per.problem
  err.mat.list[[chunk.i]] <- err.vec
}
chunk.best <- do.call(rbind, chunk.best.list)
err.mat <- do.call(rbind, err.mat.list)
err.vec <- colSums(err.mat)
train.errors <-
  data.frame(bases.per.problem=as.integer(names(err.vec)),
             errors=as.integer(err.vec),
             regions=sum(chunk.best$regions))
chosen.i <- max(which(err.vec == min(err.vec)))
res.str <- names(err.vec)[chosen.i]
bases.per.problem <- as.integer(res.str)
train.errors.picked <- train.errors[chosen.i, ]

chunk.best$selected.errors <- err.mat[, res.str]
chunk.best$best.bases.per.problem <- NA
for(chunk.i in seq_along(bpp.list)){
  bpp <- bpp.list[[chunk.i]]
  bpp[bpp == res.str] <-
    paste0("<b>", bpp[bpp==res.str], "</b>")
  chunk.best$best.bases.per.problem[chunk.i] <- paste0(bpp, collapse="<br />")
}

chunk.ordered <-
  chunk.best[order(chunk.best$selected.errors, decreasing=TRUE),]

resCurve <- 
ggplot()+
  geom_line(aes(bases.per.problem, errors),
            data=train.errors)+
  scale_x_log10()+
  ylab("minimum incorrect labels (train error)")+
  ggtitle(paste(res.str, "bases/problem"))+
  geom_point(aes(bases.per.problem, errors),
             data=train.errors.picked,
             pch=1)
print(train.errors)
cat("Picked the following problem size:\n")
print(train.errors.picked)

png.name <-
  file.path(chunks.dir, "figure-train-errors", "figure-train-errors.png")
png.dir <- dirname(png.name)
dir.create(png.dir, showWarnings=FALSE, recursive=TRUE)
png(png.name, width=400, h=300, units="px")
print(resCurve)
dev.off()

res.xt <- xtable(train.errors)
res.html <-
  print(res.xt, type="html",
        include.rownames=FALSE)

xt <- xtable(chunk.ordered)
html.table <-
  print(xt, type="html",
        include.rownames=FALSE,
        sanitize.text.function=identity)
html.out <-
  paste("<h1>Train error totals per problem size</h1>",
        "<table><tr>",
        '<td> <img src="figure-train-errors.png" /> </td>',
        "<td>", res.html, "</td>",
        "</tr></table>",
        "<h1>Train error details for each chunk of labels</h1>",
        html.table,
        "<h1>Train/validation error curves for selecting regularization</h1>",
        '<p><img src="figure-regularization.png" /></td>')

out.file <- file.path(chunks.dir, "figure-train-errors", "index.html")

cat(html.out, file=out.file)

problems.by.chunk <- list()
labeled.problems.by.chunk <- list()
regions.by.chunk <- list()
for(chunk.i in seq_along(problems.RData.vec)){
  problems.RData <- problems.RData.vec[[chunk.i]]
  objs <- load(problems.RData)
  if(! "step2.data.list" %in% objs){
    stop("step.data.list not found in ", problems.RData)
  }
  res.data <- step2.data.list[[res.str]]
  for(problem.i in 1:nrow(res.data$problems)){
    prob.info <- res.data$problems[problem.i, ]
    problem.name <- paste(prob.info$problem.name)
    target <- step2.error.list[[paste(res.str, problem.name)]]$problem$target
    mlist <- step2.model.list[[problem.name]]
    stopifnot(prob.info$problemStart < mlist$peaks$chromStart)
    stopifnot(mlist$peaks$chromEnd < prob.info$problemEnd)
    mlist$target <- target
    problems.by.chunk[[problems.RData]][[problem.name]] <- mlist
    if(is.numeric(target)){
      labeled.problems.by.chunk[[problems.RData]][[problem.name]] <- mlist
    }
  }
  chunk.dir <- dirname(problems.RData)
  regions.RData <- file.path(chunk.dir, "regions.RData")
  load(regions.RData)
  regions.by.chunk[[problems.RData]] <- regions
}

## Fit model to all training data. need (unlabeled) problems.by.chunk
## here since the predictions must be made on all problems in the
## validation set (not just the subset of problems for which we have
## assigned some regions and computed a target).
full.curves <- tv.curves(problems.by.chunk, regions.by.chunk)
picked.error <- best.on.validation(full.curves)

tvPlot <- 
  ggplot()+
    geom_point(aes(-log10(regularization), metric.value, color=tv),
               pch=1,
               data=picked.error)+
    geom_line(aes(-log10(regularization), metric.value, color=tv),
              data=full.curves)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(metric.name ~ validation.fold, labeller=function(var, val){
      if(var=="validation.fold"){
        paste("fold", val)
      }else{
        paste(val)
      }
    }, scales="free")
reg.png <-
  file.path(chunks.dir, "figure-train-errors", "figure-regularization.png")
png(reg.png, width=600, h=400, units="px")
print(tvPlot)
dev.off()

mean.reg <- mean(picked.error$regularization)

## labeled.problems.by.chunk needed here since
## IntervalRegressionProblems does not accept problems without
## targets.
train.list <- do.call(c, labeled.problems.by.chunk)
full.fit <-
  IntervalRegressionProblems(train.list,
                             initial.regularization=mean.reg,
                             factor.regularization=NULL,
                             verbose=0)

## get chrom size info from a bigwig, so we can generate a list of
## segmentation problems and divide them into jobs.
bigwig.file.vec <- Sys.glob(file.path(data.dir, "*", "*.bigwig"))
bigwig.file <- bigwig.file.vec[1]
chrom.ranges <- bigWigInfo(bigwig.file)
ranges.by.chrom <- split(chrom.ranges, chrom.ranges$chrom)
problems.by.chrom <- list()
for(chrom in names(ranges.by.chrom)){
  chrom.range <- ranges.by.chrom[[chrom]]
  all.chrom.problems <- with(chrom.range, {
    getProblems(chrom, chromStart, chromEnd, bases.per.problem,
                chrom.size=chromEnd)
  })
  ## Look at a bigwig file to see where the first chromStart and last
  ## chromEnd are, and then only process the problems which have some
  ## data.
  cmd <-
    sprintf("bigWigToBedGraph %s /dev/stdout -chrom=%s",
            bigwig.file, chrom)
  head.cmd <- paste(cmd, "| head -1")
  tail.cmd <- paste(cmd, "| tail -1")
  head.dt <- fread.or.null(head.cmd)
  tail.dt <- fread.or.null(tail.cmd)
  two <- rbind(head.dt, tail.dt)
  setnames(two, c("chrom", "chromStart", "chromEnd", "count"))
  first.chromStart <- two$chromStart[1]
  last.chromEnd <- two$chromEnd[2]
  problems.by.chrom[[chrom]] <- 
    all.chrom.problems[problemStart < last.chromEnd &
                         first.chromStart < problemEnd, ]
}
overlapping.problems <- do.call(rbind, problems.by.chrom)
overlapping.problems$job <- sort(rep(1:numJobs, l=nrow(overlapping.problems)))
table(overlapping.problems$job)
problems.by.job <- split(overlapping.problems, overlapping.problems$job)

trained.model.RData <- file.path(chunks.dir, "trained.model.RData")
save(train.errors, train.errors.picked,
     full.fit,
     ## for estimating test error later:
     labeled.problems.by.chunk,
     problems.by.chunk,
     regions.by.chunk,
     ## for parallelizing prediction on jobs:
     problems.by.job,
     file=trained.model.RData)

