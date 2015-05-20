library(animint)
require(data.table)
require(ggplot2)

argv <- "~/exampleData/PeakSegJoint-chunks/1"

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

if(length(argv) != 1){
  stop("usage: Step4.R path/to/PeakSegJoint-chunks/1")
}

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

chunk.dir <- argv
chunks.dir <- dirname(chunk.dir)
trained.model.RData <- file.path(chunks.dir, "trained.model.RData")
problems.RData <- file.path(chunk.dir, "problems.RData")
regions.RData <- file.path(chunk.dir, "regions.RData")

robjs <- load(regions.RData)
regions$region.i <- 1:nrow(regions)
chrom <- paste(regions$chrom[1])

## Filter counts by limits (area around regions).
limits <- regions[, .(chromStart=min(chromStart),
                      chromEnd=max(chromEnd))]
limits[, bases := chromEnd - chromStart ]
limits[, expand := as.integer(bases/10) ]
limits[, min := chromStart - expand]
limits[, max := chromEnd + expand]
lim.vec <- with(limits, c(min, max))/1e3

mobjs <- load(trained.model.RData)
pobjs <- load(problems.RData)

## Some fake data to demonstrate how to sub-sample a coverage
## data.frame using approx.
fake.segs <- data.frame(chromStart=c(0, 1, 3),
                   chromEnd=c(1, 3, 10),
                   coverage=1:3)
x <- 1:10
y <- with(fake.segs, {
  approx(chromStart+1, coverage, x,
         method="constant", rule=2)
})$y
fake.points <- data.frame(x,y)
## ggplot()+
##   geom_segment(aes(chromStart+0.5, coverage,
##                    xend=chromEnd+0.5, yend=coverage),
##                data=fake.segs)+
##   geom_point(aes(x, y), data=fake.points)

## Read and subsample count data to plot width.
width.pixels <- 1000
counts.RData.vec <- Sys.glob(file.path(chunk.dir, "*", "*.RData"))
counts.by.sample <- list()
for(counts.RData.path in counts.RData.vec){
  objs <- load(counts.RData.path)
  sample.id <- sub(".RData$", "", basename(counts.RData.path))
  a.data <- with(counts, {
    approx(chromStart+1, count,
           as.integer(seq(limits$min, limits$max, l=width.pixels)),
           method="constant", rule=2)
  })
  counts.by.sample[[sample.id]] <- with(a.data, {
    data.table(sample.id,
               base=x,
               count=y)
  })
}
some.counts <- do.call(rbind, counts.by.sample)

step2.problems.by.res <- list()
prob.labels.by.res <- list()
for(res.str in names(step2.data.list)){
  res.data <- step2.data.list[[res.str]]
  step2.problems.by.res[[res.str]] <- 
    data.table(sample.id="problems", res.data$problems)
  bases.vec <- with(res.data$problems, problemEnd-problemStart)
  prob.labels.by.res[[res.str]] <-
    data.table(sample.id="problems",
               bases.per.problem=as.integer(res.str),
               mean.bases=as.integer(mean(bases.vec)),
               problems=nrow(res.data$problems))
}
step2.problems <- do.call(rbind, step2.problems.by.res)
prob.labels <- do.call(rbind, prob.labels.by.res)
prob.labels$problem.i <- max(prob.labels$problems)
prob.labels$chromStart <- limits$chromStart

facet.rows <- length(counts.by.sample)+1
dvec <- diff(log(res.error$bases.per.problem))
dval <- exp(mean(dvec))
dval2 <- (dval-1)/2 + 1
viz <-
  list(coverage=ggplot()+
         geom_text(aes(chromStart/1e3, problem.i,
                       showSelected=bases.per.problem,
                       label=sprintf("%d problems mean size %.1f kb",
                         problems, mean.bases/1e3)),
                   data=prob.labels,
                   hjust=0)+
         geom_segment(aes(problemStart/1e3, problem.i,
                          showSelected=bases.per.problem,
                          xend=problemEnd/1e3, yend=problem.i),
                      data=step2.problems)+
         scale_y_continuous("aligned read coverage",
                            breaks=function(limits){
                              floor(limits[2])
                            })+
         scale_linetype_manual("error type",
                               limits=c("correct", 
                                 "false negative",
                                 "false positive"
                                        ),
                               values=c(correct=0,
                                 "false negative"=3,
                                 "false positive"=1))+
         scale_x_continuous(paste("position on",
                                  chrom,
                                  "(kilo bases = kb)"))+
         coord_cartesian(xlim=lim.vec)+
         geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                           fill=annotation),
                       alpha=0.5,
                       color="grey",
                       data=regions)+
         scale_fill_manual(values=ann.colors)+
         theme_bw()+
         theme_animint(width=width.pixels, height=facet.rows*100)+
         theme(panel.margin=grid::unit(0, "cm"))+
         facet_grid(sample.id ~ ., labeller=function(var, val){
           sub("McGill0", "", sub(" ", "\n", val))
         }, scales="free")+
         geom_line(aes(base/1e3, count),
                   data=some.counts,
                   color="grey50"),

       resError=ggplot()+
           geom_tallrect(aes(xmin=bases.per.problem/dval2,
                             xmax=bases.per.problem*dval2,
                             clickSelects=bases.per.problem),
                         alpha=0.5,
                         data=res.error)+
           scale_x_log10()+
           geom_line(aes(bases.per.problem, errors/regions*100,
                         color=chunks, size=chunks),
                     data=data.frame(res.error, chunks="this"))+
           geom_line(aes(bases.per.problem, errors/regions*100,
                         color=chunks, size=chunks),
                     data=data.frame(train.errors, chunks="all")),
                     
       title=chunk.dir,

       first=list(bases.per.problem=train.errors.picked$bases.per.problem))

animint.dir <- file.path(chunk.dir, "figure-train-errors")
animint2dir(viz, animint.dir)

