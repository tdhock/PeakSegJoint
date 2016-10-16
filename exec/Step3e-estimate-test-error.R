library(PeakSegJoint)

model.path <-
  system.file("exampleData",
              "PeakSegJoint-chunks",
              "trained.model.RData",
              package="PeakSegJoint")
model.path <- "~/exampleData/PeakSegJoint-chunks/trained.model.RData"
model.path <- "~/projects/H3K27ac_TDH/PeakSegJoint-chunks/trained.model.RData"
argv <- c(model.path, "1", "1")

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

ppn <- PPN.cores()
if(!is.na(ppn))options(mc.cores=ppn/2)

if(length(argv) != 3){
  stop("usage: Step3e.R path/to/PeakSegJoint-chunks/trained.model.RData chunk.order.seed test.fold")
}

trained.model.RData <- normalizePath(argv[1], mustWork=TRUE)
chunk.order.seed <- as.integer(argv[2])
test.fold <- as.integer(argv[3])
chunks.dir <- dirname(trained.model.RData)
data.dir <- dirname(chunks.dir)

model.objs <- load(trained.model.RData)

## Divide chunks into train+validation/test.
map.RData <- file.path(data.dir, "chunk.file.map.RData")
load(map.RData)
chunks.by.file <- split(chunk.file.map, chunk.file.map$labels.file)
all.chunk.names <- names(problems.by.chunk)
names(all.chunk.names) <- basename(dirname(all.chunk.names))
cv.chunk.names <- all.chunk.names[paste(chunk.file.map$chunk.id)]

labeled.problems.by.chunk <- list()
for(chunk.name in names(problems.by.chunk)){
  by.problem <- problems.by.chunk[[chunk.name]]
  has.target <- sapply(by.problem, function(L)is.numeric(L$target))
  labeled.problems.by.chunk[[chunk.name]] <- by.problem[has.target]
}

## For each test fold, hold it out and train a sequence of models with
## increasingly more data, and compute test error of each.
test.metrics.curve.list <- list()

is.test <- outer.fold.id == test.fold
set.seed(chunk.order.seed)
sets <-
  list(test=cv.chunk.names[is.test],
       train.validation=sample(cv.chunk.names[!is.test]))
test.regions <- regions.by.chunk[sets$test]
test.problems <- problems.by.chunk[sets$test]

if(length(sets$train.validation) < 2){
  print(sets$train.validation)
  stop("need at least 2 train chunks, please add more labels")
}
train.chunks.vec <- 2:length(sets$train.validation)
for(train.chunks in train.chunks.vec){
  print(system.time({
    cat("estimating test error:",
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
                                 factor.regularization=NULL,
                                 verbose=0)
    test.results <- error.metrics(test.problems, test.regions, tv.fit)
    test.results$metrics$test.fold <- test.fold
    test.results$metrics$chunk.order.seed <- chunk.order.seed
    test.results$metrics$train.chunks <- train.chunks
    stopifnot(test.results$metrics["incorrect.regions", "possible"] ==
                sum(sapply(test.regions, nrow)))
    test.metrics.curve.list[[paste(train.chunks)]] <- test.results$metrics
  }))
}
test.metrics.curve <- do.call(rbind, test.metrics.curve.list)
rownames(test.metrics.curve) <- NULL

test.out.dir <- file.path(chunks.dir, "figure-test-errors")
seed.RData <-
  file.path(test.out.dir,
            paste0("seed", chunk.order.seed, "fold", test.fold, ".RData"))

save(test.results,
     test.metrics.curve,
     file=seed.RData)

## ggplot()+
##   ggtitle(paste("test error for one",
##                 "random ordering of the labeled train chunks"))+
##   geom_text(aes(train.chunks, metric.value/possible*100,
##                 label=sprintf("%.1f%%", metric.value/possible*100)),
##             data=test.results$metrics,
##             vjust=-1)+
##   ylab("percent incorrect (test error)")+
##   scale_x_continuous("number of labeled chunks in train set",
##                      breaks=function(lim.vec){
##                        ##print(lim.vec)
##                        ceiling(lim.vec[1]):floor(lim.vec[2])
##                      })+
##   geom_line(aes(train.chunks, metric.value/possible*100,
##                 group=chunk.order.seed),
##             data=test.metrics.curve)+
##   geom_point(aes(train.chunks, metric.value/possible*100,
##                  group=chunk.order.seed),
##              data=test.metrics.curve)+
##   theme_bw()+
##   theme(panel.margin=grid::unit(0, "cm"))+
##   facet_grid(metric.name ~ test.fold, scales="free")

