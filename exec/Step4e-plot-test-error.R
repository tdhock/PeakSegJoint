library(PeakSegJoint)
library(ggplot2)
library(xtable)
options(xtable.print.results=FALSE)

argv <- system.file("exampleData", "PeakSegJoint-chunks",
                    mustWork=TRUE,
                    package="PeakSegJoint")

argv <- commandArgs(trailingOnly=TRUE)

if(length(argv) != 1){
  stop("Usage: Step4e.R path/to/PeakSegJoint-chunks")
}

chunks.dir <- normalizePath(argv[1])
data.dir <- dirname(chunks.dir)

regions.RData.vec <- Sys.glob(file.path(chunks.dir, "*", "regions.RData"))
regions.list <- list()
for(regions.RData in regions.RData.vec){
  load(regions.RData)
  chunk.id.str <- basename(dirname(regions.RData))
  regions.list[[chunk.id.str]] <- regions
}

map.objs <- load(file.path(data.dir, "chunk.file.map.RData"))

test.out.dir <- file.path(chunks.dir, "figure-test-errors")
RData.file.vec <- Sys.glob(file.path(test.out.dir, "seed*.RData"))
test.regions.list <- list()
test.peaks.list <- list()
test.folds.list <- list()
test.metrics.curves.list <- list()
test.metrics.most.list <- list()
for(RData.file in RData.file.vec){
  objs <- load(RData.file)
  test.fold <- test.metrics.curve$test.fold[1]
  test.metrics.most.list[[paste(test.fold)]] <- test.results$metrics
  test.metrics.curves.list[[RData.file]] <- test.metrics.curve
  test.peaks.list[names(test.results$peaks)] <- test.results$peaks
  for(chunk.name in names(test.results$error.regions)){
    error.regions <- test.results$error.regions[[chunk.name]]
    chunk.id.str <- basename(dirname(chunk.name))
    orig.regions <- regions.list[[chunk.id.str]]
    stopifnot(nrow(orig.regions) == nrow(error.regions))
    ##stopifnot(nrow(error.regions) == test.results$)
    test.regions.list[[chunk.name]] <- error.regions
    test.folds.list[[chunk.name]] <- with(error.regions, {
      data.table(chrom=chrom[1],
                 test.fold,
                 position=(min(chromStart)+max(chromEnd))/2)
    })
  }
}
test.folds <- do.call(rbind, test.folds.list)
test.metrics.most <- do.call(rbind, test.metrics.most.list)
test.metrics.curves <- do.call(rbind, test.metrics.curves.list)
rownames(test.metrics.curves) <- NULL

## get chrom ranges for plot.
bigwig.file.vec <- Sys.glob(file.path(data.dir, "*", "*.bigwig"))
bigwig.file <- bigwig.file.vec[1]
chrom.ranges <- bigWigInfo(bigwig.file)

chunk.counts <- table(test.folds$test.fold)
ggfolds <- ggplot()+
  ggtitle("Distribution of folds across chromosomes")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(test.fold ~ ., labeller=function(var, val){
    n.chunks <- chunk.counts[paste(val)]
    paste0("test fold ", val, "\n", n.chunks, " chunks")
  }, scales="free", space="free")+
  geom_segment(aes(chromStart/1e6, chrom,
                   xend=chromEnd/1e6, yend=chrom),
               data=chrom.ranges)+
  geom_point(aes(position/1e6, chrom),
             pch=1,
             data=test.folds)+
  xlab("position on chromosome (mega bases = Mb)")
folds.png <-
  file.path(test.out.dir, "figure-folds.png")
png(folds.png, width=1000, h=600, units="px")
print(ggfolds)
dev.off()

gg.test <- 
  ggplot()+
    ggtitle(paste("test error for one",
                  "random ordering of the labeled train chunks"))+
    geom_text(aes(train.chunks, metric.value/possible*100,
                  label=sprintf("%.1f%%", metric.value/possible*100)),
              data=test.metrics.most,
              vjust=-1)+
    ylab("percent incorrect (test error)")+
    scale_x_continuous("number of labeled chunks in train set",
                       breaks=function(lim.vec){
                         ##print(lim.vec)
                         ceiling(lim.vec[1]):floor(lim.vec[2])
                       })+
    geom_line(aes(train.chunks, metric.value/possible*100,
                  group=chunk.order.seed),
              data=test.metrics.curves)+
    geom_point(aes(train.chunks, metric.value/possible*100,
                  group=chunk.order.seed),
              data=test.metrics.curves)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(metric.name ~ test.fold, labeller=function(var, val){
      if(var=="test.fold"){
        paste("test fold", val)
      }else{
        paste(val)
      }
    }, scales="free")

test.png <-
  file.path(test.out.dir, "figure-test-error-decreases.png")
png(test.png, width=1000, h=600, units="px")
print(gg.test)
dev.off()

most.list <- split(test.metrics.most, test.metrics.most$test.fold)
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
        '<img src="figure-folds.png" alt="folds" />',
        '<br />',
        '<img src="figure-test-error-decreases.png" alt="error" />',
        "<h1>Test error details for each chunk of labels</h1>",
        test.html.table)
cat(test.html.out, file=test.index.html)

save(test.metrics.curves,
     test.regions.list,
     test.error.summary,
     file=file.path(test.out.dir, "test.metrics.curves.RData"))

totals <-
  rbind(regions=sum(sapply(regions.list, nrow)),
        test.regions=sum(sapply(test.regions.list, nrow)),
        metrics=test.error.summary["incorrect.regions", "possible"])
print(totals)

