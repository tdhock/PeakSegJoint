library(ggplot2)
library(xtable)
library(PeakSegJoint)
options(xtable.print.results=FALSE)

argv <-
  system.file("exampleData",
              "PeakSegJoint-chunks",
              "trained.model.RData",
              package="PeakSegJoint")

argv <- "~/exampleData/PeakSegJoint-chunks/trained.model.RData"
argv <- "~/projects/H3K27ac_TDH/PeakSegJoint-chunks/trained.model.RData"

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

PPN.cores()

if(length(argv) != 1){
  stop("usage: Step3e.R path/to/PeakSegJoint-chunks/trained.model.RData")
}

trained.model.RData <- normalizePath(argv)
chunks.dir <- dirname(trained.model.RData)
data.dir <- dirname(chunks.dir)

load(trained.model.RData)

## Divide chunks into train+validation/test.
map.RData <- file.path(data.dir, "chunk.file.map.RData")
load(map.RData)
chunks.by.file <- split(chunk.file.map, chunk.file.map$labels.file)
all.chunk.names <- names(problems.by.chunk)
names(all.chunk.names) <- basename(dirname(all.chunk.names))
if(length(chunks.by.file) == 1){
  outer.folds <- 4
  if(length(all.chunk.names) < outer.folds){
    outer.folds <- length(all.chunk.names)
  }
  outer.fold.id <- sample(rep(1:outer.folds, l=length(all.chunk.names)))
  names(outer.fold.id) <- names(all.chunk.names)
  fold.msg <- "randomly selected folds"
}else{
  outer.folds <- length(chunks.by.file)
  outer.fold.id <- rep(NA, l=length(all.chunk.names))
  names(outer.fold.id) <- names(all.chunk.names)
  for(file.i in seq_along(chunks.by.file)){
    file.chunks <- chunks.by.file[[file.i]]
    outer.fold.id[paste(file.chunks$chunk.id)] <- file.i
  }
  fold.msg <- "one fold for each labels file"
}
test.error.msg <-
  paste0(outer.folds, " fold cross-validation (", fold.msg, ").")

## For each test fold, hold it out and train a sequence of models with
## increasingly more data, and compute test error of each.
test.error.list <- list()
test.peaks.list <- list()
test.regions.list <- list()
for(test.fold in 1:outer.folds){
  is.test <- outer.fold.id == test.fold
  n.chunk.order.seeds <- 2
  for(chunk.order.seed in 1:n.chunk.order.seeds){
    set.seed(chunk.order.seed)
    sets <- list(train.validation=sample(all.chunk.names[!is.test]),
                 test=all.chunk.names[is.test])
    if(length(sets$train.validation) < 2){
      print(sets$train.validation)
      stop("need at least 2 train chunks, please add more labels")
    }
    train.chunks.vec <- 2:length(sets$train.validation)
    for(train.chunks in train.chunks.vec){
      cat("estimating test error:",
          test.fold, "/", outer.folds, "folds,",
          chunk.order.seed, "/", n.chunk.order.seeds, "seeds,",
          train.chunks, "/", length(sets$train.validation), "chunks.\n")
      some.train.validation <- sets$train.validation[1:train.chunks]
      some.problems <- problems.by.chunk[some.train.validation]
      some.regions <- regions.by.chunk[some.train.validation]
      full.curves <- tv.curves(some.problems, some.regions)
      picked.error <- best.on.validation(full.curves)
      mean.reg <- mean(picked.error$regularization)
      tv.list <- do.call(c, labeled.problems.by.chunk[some.train.validation])
      tv.fit <-
        IntervalRegressionProblems(tv.list,
                                   initial.regularization=mean.reg,
                                   factor.regularization=10000,
                                   verbose=0)
      test.results <- error.metrics(some.problems, some.regions, tv.fit)
      test.regions.list[names(test.results$error.regions)] <-
        test.results$error.regions
      test.peaks.list[names(test.results$peaks)] <- test.results$peaks
      test.metrics <-
        subset(test.results$metrics, regularization == regularization[1])
      rownames(test.metrics) <- NULL
      test.error.list[[paste(train.chunks, chunk.order.seed, test.fold)]] <-
        data.frame(test.fold, chunk.order.seed, train.chunks, test.metrics)
    }
  }
}
test.error.curves <- do.call(rbind, test.error.list)
most.chunks <-
  subset(test.error.curves,
         chunk.order.seed==1 & train.chunks==max(train.chunks))

gg.test <- 
  ggplot()+
    ggtitle(paste("test error for",
                  n.chunk.order.seeds,
                  "random orderings of the labeled train chunks"))+
    geom_text(aes(train.chunks, metric.value/possible*100,
                  label=sprintf("%.1f%%", metric.value/possible*100)),
              data=most.chunks,
              vjust=-1)+
    ylab("percent incorrect (test error)")+
    scale_x_continuous("number of labeled chunks in train set",
                       breaks=function(lim.vec){
                         ##print(lim.vec)
                         ceiling(lim.vec[1]):floor(lim.vec[2])
                       })+
    geom_line(aes(train.chunks, metric.value/possible*100,
                  group=chunk.order.seed),
              data=test.error.curves)+
    geom_point(aes(train.chunks, metric.value/possible*100,
                  group=chunk.order.seed),
              data=test.error.curves)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(metric.name ~ test.fold, labeller=function(var, val){
      if(var=="test.fold"){
        paste("test fold", val)
      }else{
        paste(val)
      }
    }, scales="free")
test.out.dir <- file.path(chunks.dir, "figure-test-errors")
dir.create(test.out.dir, showWarnings=FALSE)
test.png <-
  file.path(test.out.dir, "figure-test-error-decreases.png")
png(test.png, width=1000, h=600, units="px")
print(gg.test)
dev.off()

most.list <- split(most.chunks, most.chunks$test.fold)
incorrect <- as.integer(rowSums(sapply(most.list, "[[", "metric.value")))
possible <- as.integer(rowSums(sapply(most.list, "[[", "possible")))
percent.incorrect <- incorrect / possible * 100
test.error.summary <- 
  data.frame(row.names=most.list[[1]]$metric.name,
             incorrect, possible, percent.incorrect)

test.index.html <- file.path(test.out.dir, "index.html")
test.row.list <- list()
for(problems.RData in names(test.regions.list)){
  chunk.dir <- dirname(problems.RData)
  chunk.id <- basename(chunk.dir)
  test.error.figure <-
    sprintf('<img src="%s.png" alt="chunk%s" />', chunk.id, chunk.id)
  error.regions <- test.regions.list[[problems.RData]]
  chunk.peaks <- test.peaks.list[[problems.RData]]
  out.RData <- file.path(test.out.dir, paste0(chunk.id, ".RData"))
  save(error.regions, chunk.peaks, file=out.RData)
  test.row.list[[chunk.id]] <- with(error.regions, {
    data.frame(test.error.figure,
               errors=sum(fp+fn),
               regions=length(fp),
               fp=sum(fp),
               possible.fp=sum(possible.fp),
               fn=sum(fn),
               possible.fn=sum(possible.tp),
               test.fold=outer.fold.id[[chunk.id]])
  })
}
test.row.df <- do.call(rbind, test.row.list)
test.row.df <- test.row.df[order(test.row.df$errors, decreasing=TRUE), ]
test.xt <- xtable(test.row.df)
test.html.table <-
  print(test.xt, type="html",
        include.rownames=FALSE,
        sanitize.text.function=identity)
summary.xt <- xtable(test.error.summary)
summary.html <- print(summary.xt, type="html",
                      include.rownames=TRUE)
test.html.out <-
  paste("<h1>Test error summary</h1>",
        summary.html,
        "<p>Targets counts examples (genomic regions to segment)",
        "in the interval regression problem.</p>",
        "<p>Regions, false postives, and false negatives",
        "count labels (peakStart, peakEnd, peaks, noPeaks).</p>",
        "<p>Test error was estimated using",
        test.error.msg,
        "</p>",
        '<img src="figure-test-error-decreases.png" />',
        "<h1>Test error details for each chunk of labels</h1>",
        test.html.table)
cat(test.html.out, file=test.index.html)

stopifnot(test.error.summary["incorrect.regions", "possible"] ==
            sum(sapply(regions.by.chunk, nrow)))

save(test.error.curves,
     file=file.path(test.out.dir, "test.error.curves.RData"))
