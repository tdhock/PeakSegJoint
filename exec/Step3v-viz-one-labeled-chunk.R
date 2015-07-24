library(animint)
library(PeakSegJoint)

argv <- "~/exampleData/PeakSegJoint-chunks/3"

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

if(length(argv) != 1){
  stop("usage: Step3v.R path/to/PeakSegJoint-chunks/1")
}

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

chunk.dir <- argv
chunks.dir <- dirname(chunk.dir)
data.dir <- dirname(chunks.dir)
trained.model.RData <- file.path(chunks.dir, "trained.model.RData")
problems.RData <- file.path(chunk.dir, "problems.RData")
regions.RData <- file.path(chunk.dir, "regions.RData")

mobjs <- load(trained.model.RData)
pobjs <- load(problems.RData)
robjs <- load(regions.RData)
regions$region.i <- 1:nrow(regions)
chrom <- paste(regions$chrom[1])

width.pixels <- 1500
bigwig.file.vec <- Sys.glob(file.path(data.dir, "*", "*.bigwig"))
limits.by.sample <- list()
for(bigwig.file in bigwig.file.vec){
  sample.counts <-
    readBigWig(bigwig.file, chunk$chrom, chunk$chunkStart, chunk$chunkEnd)
  sample.id <- sub("[.]bigwig$", "", basename(bigwig.file))
  limits.by.sample[[sample.id]] <- with(sample.counts, {
    data.table(sample.id,
               min.chromStart=min(chromStart),
               max.chromEnd=max(chromEnd))
  })
}

region.limits <- regions[, .(chromStart=min(chromStart),
                             chromEnd=max(chromEnd))]
region.limits[, bases := chromEnd - chromStart ]
region.limits[, expand := as.integer(bases/10) ]
region.limits[, min := chromStart - expand]
region.limits[, max := chromEnd + expand]
limits.by.sample[["regions"]] <- with(region.limits, {
  data.table(sample.id="regions",
             min.chromStart=min,
             max.chromEnd=max)
})
sample.limits <- do.call(rbind, limits.by.sample)

## We can't simply use chunkStart/chunkEnd, because that is much too
## zoomed out (the big problem sizes start/end much outside the
## first/last annotated region labels.
limits <- sample.limits[, .(min=max(min.chromStart),
                            max=min(max.chromEnd))]
lim.vec <- with(limits, c(min, max))/1e3

counts.by.sample <- list()
for(bigwig.file in bigwig.file.vec){
  sample.counts <-
    readBigWig(bigwig.file, chunk$chrom, chunk$chunkStart, chunk$chunkEnd)
  sample.id <- sub("[.]bigwig$", "", basename(bigwig.file))
  a.data <- with(sample.counts, {
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
modelSelection.by.problem <- list()
peaks.by.problem <- list()
first.selection.list <-
  list(bases.per.problem=train.errors.picked$bases.per.problem)
regions.by.problem <- list()
for(res.str in names(step2.data.list)){
  res.data <- step2.data.list[[res.str]]
  step2.problems.by.res[[res.str]] <- 
    data.table(sample.id="problems", res.data$problems)
  bases.vec <- with(res.data$problems, problemEnd-problemStart)
  bases.per.problem <- as.integer(res.str)
  prob.labels.by.res[[res.str]] <-
    data.table(sample.id="problems",
               bases.per.problem,
               mean.bases=as.integer(mean(bases.vec)),
               problems=nrow(res.data$problems))
  for(problem.i in 1:nrow(res.data$problems)){
    problem <- res.data$problems[problem.i, ]
    problem.name <- paste(problem$problem.name)
    problem.dot <- paste0(gsub("[:-]", ".", problem.name), "peaks")
    error <- step2.error.list[[paste(res.str, problem.name)]]
    model <- step2.model.list[[problem.name]]
    if(!is.null(model)){
      if(is.data.frame(model$peaks)){
        peaks.by.problem[[paste(res.str, problem.dot)]] <-
          data.table(problem, model$peaks)
      }
      if(is.null(error$peaks)){
        first.selection.list[[problem.dot]] <- 0
      }else{
        ms <- 
        first.selection.list[[problem.dot]] <- error$peaks$peaks[1]
      }
      
      ms <- if(is.list(error$problem$error.regions)){
        regions.by.peaks <- list()
        for(peaks.str in names(error$problem$error.regions)){
          regions.df <- error$problem$error.regions[[peaks.str]]
          regions.df$peaks <- as.integer(peaks.str)
          regions.by.peaks[[peaks.str]] <-
            data.table(problem, regions.df)
        }
        regions.by.problem[[paste(res.str, problem.dot)]] <-
          do.call(rbind, regions.by.peaks)
        error$problem$modelSelection
      }else{
        data.frame(model$modelSelection,
                   errors=NA)
      }      
      modelSelection.by.problem[[paste(res.str, problem.dot)]] <-
        data.table(problem, ms)
    }#if(is.null(model
  }#for(problem.i
}
step2.problems <- do.call(rbind, step2.problems.by.res)
prob.labels <- do.call(rbind, prob.labels.by.res)
prob.labels$problem.i <- max(prob.labels$problems)
prob.labels$chromStart <- region.limits$chromStart

## prob.regions are the black segments that show which regions are
## mapped to which segmentation problems.
all.regions <- do.call(rbind, regions.by.problem)
prob.regions.names <-
  c("bases.per.problem", "problem.i", "problem.name",
    "chromStart", "chromEnd")
prob.regions <- unique(data.frame(all.regions)[, prob.regions.names])
prob.regions$sample.id <- "problems"

all.modelSelection <- do.call(rbind, modelSelection.by.problem)
modelSelection.errors <- all.modelSelection[!is.na(errors), ]
penalty.range <-
  with(all.modelSelection, c(min(max.log.lambda), max(min.log.lambda)))
penalty.mid <- mean(penalty.range)

facet.rows <- length(counts.by.sample)+1
dvec <- diff(log(res.error$bases.per.problem))
dval <- exp(mean(dvec))
dval2 <- (dval-1)/2 + 1
res.error$min.bases.per.problem <- res.error$bases.per.problem/dval2
res.error$max.bases.per.problem <- res.error$bases.per.problem*dval2

modelSelection.labels <-
  unique(all.modelSelection[, {
    list(problem.name=problem.name,
         bases.per.problem=bases.per.problem,
         problemStart=problemStart,
         problemEnd=problemEnd,
         min.log.lambda=penalty.mid,
         peaks=max(peaks)+0.5)
  }])

sample.peaks <- do.call(rbind, peaks.by.problem)
prob.peaks.names <-
  c("bases.per.problem", "problem.i", "problem.name", "peaks",
    "chromStart", "chromEnd")
problem.peaks <- unique(data.frame(sample.peaks)[, prob.peaks.names])
problem.peaks$sample.id <- "problems"

max.height.pixels <- 1000
ideal.pixels.per.facet <- 100
ideal.height <- facet.rows * ideal.pixels.per.facet
height.pixels <- as.integer(if(max.height.pixels < ideal.height){
  max.height.pixels
}else{
  ideal.height
})

peakvar <- function(position){
  paste0(gsub("[-:]", ".", position), "peaks")
}
cat("constructing data viz with .variable .value\n")
print(system.time({
  viz <-
    list(title=paste("PeakSegJoint model train errors for",
           chunk.dir),

         resError=ggplot()+
           ggtitle("select problem size")+
           ylab("minimum percent incorrect regions")+
           geom_tallrect(aes(xmin=min.bases.per.problem,
                             xmax=max.bases.per.problem,
                             clickSelects=bases.per.problem),
                         alpha=0.5,
                         data=res.error)+
           theme_animint(height=300)+
           scale_x_log10()+
           geom_line(aes(bases.per.problem, errors/regions*100,
                         color=chunks, size=chunks),
                     data=data.frame(res.error, chunks="this"))+
           geom_line(aes(bases.per.problem, errors/regions*100,
                         color=chunks, size=chunks),
                     data=data.frame(train.errors, chunks="all")),

         modelSelection=ggplot()+
           ggtitle("select number of samples with 1 peak")+
           theme_animint(height=300)+
           geom_segment(aes(min.log.lambda, peaks,
                            xend=max.log.lambda, yend=peaks,
                            showSelected=problem.name,
                            showSelected2=bases.per.problem),
                        data=data.frame(all.modelSelection, what="peaks"),
                        size=5)+
           geom_text(aes(min.log.lambda, peaks,
                         showSelected=problem.name,
                         showSelected2=bases.per.problem,
                         label=sprintf("%.1f kb in problem %s",
                           (problemEnd-problemStart)/1e3, problem.name)),
                     data=data.frame(modelSelection.labels, what="peaks"))+
           geom_segment(aes(min.log.lambda, as.integer(errors),
                            xend=max.log.lambda, yend=as.integer(errors),
                            showSelected=problem.name,
                            showSelected2=bases.per.problem),
                        data=data.frame(modelSelection.errors, what="errors"),
                        size=5)+
           ylab("")+
           xlab("penalty log.lambda")+
           geom_tallrect(aes(xmin=min.log.lambda, 
                             xmax=max.log.lambda, 
                             clickSelects.variable=
                               peakvar(problem.name),
                             clickSelects.value=peaks,
                             showSelected=problem.name,
                             showSelected2=bases.per.problem),
                         data=all.modelSelection, alpha=0.5)+
           facet_grid(what ~ ., scales="free"),
         
         coverage=ggplot()+
           ggtitle("select problem")+
           geom_segment(aes(chromStart/1e3, problem.i,
                            xend=chromEnd/1e3, yend=problem.i,
                            showSelected=bases.per.problem,
                            clickSelects=problem.name),
                        data=prob.regions)+
           geom_text(aes(chromStart/1e3, problem.i,
                         showSelected=bases.per.problem,
                         label=sprintf("%d problems mean size %.1f kb",
                           problems, mean.bases/1e3)),
                     data=prob.labels,
                     hjust=0)+
           geom_segment(aes(problemStart/1e3, problem.i,
                            showSelected=bases.per.problem,
                            clickSelects=problem.name,
                            xend=problemEnd/1e3, yend=problem.i),
                        size=5,
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
           scale_x_continuous(paste("position on", chrom,
                                    "(kilo bases = kb)"))+
           coord_cartesian(xlim=lim.vec)+
           geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                             fill=annotation),
                         alpha=0.5,
                         color="grey",
                         data=regions)+
           scale_fill_manual(values=ann.colors)+
           theme_bw()+
           theme_animint(width=width.pixels, height=height.pixels)+
           theme(panel.margin=grid::unit(0, "cm"))+
           facet_grid(sample.id ~ ., labeller=function(var, val){
             sub("McGill0", "", sub(" ", "\n", val))
           }, scales="free")+
           geom_line(aes(base/1e3, count),
                     data=some.counts,
                     color="grey50")+
           geom_tallrect(aes(xmin=chromStart/1e3,
                             xmax=chromEnd/1e3,
                             linetype=status,
                             showSelected.value=peaks,
                             showSelected.variable=peakvar(problem.name),
                             showSelected2=bases.per.problem),
                         data=all.regions,
                         fill=NA,
                         color="black")+
           geom_segment(aes(chromStart/1e3, 0,
                            xend=chromEnd/1e3, yend=0,
                            clickSelects=problem.name,
                            showSelected.variable=peakvar(problem.name),
                            showSelected.value=peaks,
                            showSelected2=bases.per.problem),
                        data=sample.peaks, size=7, color="deepskyblue")+
           geom_point(aes(chromStart/1e3, 0,
                          clickSelects=problem.name,
                          showSelected.variable=peakvar(problem.name),
                          showSelected.value=peaks,
                          showSelected2=bases.per.problem),
                      data=sample.peaks, size=3, color="deepskyblue")+
           geom_segment(aes(chromStart/1e3, problem.i,
                            xend=chromEnd/1e3, yend=problem.i,
                            clickSelects=problem.name,
                            showSelected.variable=peakvar(problem.name),
                            showSelected.value=peaks,
                            showSelected2=bases.per.problem),
                        data=problem.peaks, size=7, color="deepskyblue"),

         first=first.selection.list)

  ## For every problem there is a selector (called problem.name) for
  ## the number of peaks in that problem. There is also a selection
  ## variable for every unique value of clickSelects.variable and
  ## showSelected.variable.
}))

animint.dir <- file.path(chunk.dir, "figure-train-errors")

cat("compiling data viz\n")
print(system.time({
  animint2dir(viz, out.dir=animint.dir)
}))

first.peaks.all <-
  sample.peaks[bases.per.problem == first.selection.list$bases.per.problem, ]
first.peaks.all[, selector.name := peakvar(problem.name) ]
setkey(first.peaks.all, selector.name, peaks)
some.first <- first.selection.list[unique(first.peaks.all$selector.name)]
first.dt <- data.table(selector.name=names(some.first),
                       peaks=unlist(some.first))
first.peaks <- first.peaks.all[first.dt, ]
peaks.not.na <- first.peaks[!is.na(sample.id),]
all.regions[, selector.name := peakvar(problem.name)]
setkey(all.regions, selector.name, peaks, bases.per.problem)
first.dt$bases.per.problem <- first.selection.list$bases.per.problem
first.regions <- all.regions[first.dt, ]
regions.not.na <- first.regions[!is.na(sample.id), ]

## This test was originally designed to make sure that the number of
## errors that we show is the same as the number of regions, under the
## assumption that each region is assigned to exactly 1 segmentation
## problem. However in some rare cases some regions are assigned to a
## segmentation problem for which no peaks were feasible. In this case
## we will not show the error regions, since anyways it does not help
## us train the model.

##stopifnot(nrow(regions.not.na) == nrow(regions))

train.fig <- 
  ggplot()+
    ggtitle(animint.dir)+
    scale_y_continuous("aligned read coverage",
                       breaks=function(limits){
                         floor(limits[2])
                       })+
    scale_x_continuous(paste("position on", chrom,
                             "(kilo bases = kb)"))+
    coord_cartesian(xlim=lim.vec)+
    geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                      fill=annotation),
                  alpha=0.5,
                  color="grey",
                  data=regions)+
           scale_linetype_manual("error type",
                                 values=c(correct=0,
                                   "false negative"=3,
                                   "false positive"=1))+
    scale_fill_manual(values=ann.colors)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., labeller=function(var, val){
      sub("McGill0", "", sub(" ", "\n", val))
    }, scales="free")+
    geom_line(aes(base/1e3, count),
              data=some.counts,
              color="grey50")+
    geom_tallrect(aes(xmin=chromStart/1e3,
                      xmax=chromEnd/1e3,
                      linetype=status,
                      showSelected.value=peaks,
                      showSelected.variable=peakvar(problem.name),
                      showSelected2=bases.per.problem),
                  data=regions.not.na,
                  fill=NA,
                  size=1.5,
                  color="black")+
    geom_segment(aes(chromStart/1e3, 0,
                     xend=chromEnd/1e3, yend=0,
                     clickSelects=problem.name,
                     showSelected.variable=peakvar(problem.name),
                     showSelected.value=peaks,
                     showSelected2=bases.per.problem),
                 data=peaks.not.na, size=4, color="deepskyblue")+
    geom_point(aes(chromStart/1e3, 0,
                   clickSelects=problem.name,
                   showSelected.variable=peakvar(problem.name),
                   showSelected.value=peaks,
                   showSelected2=bases.per.problem),
               data=peaks.not.na, size=4, color="deepskyblue")

png.name <- paste0(animint.dir, ".png")
png(png.name, width=1000, h=(facet.rows+1)*30, units="px")
print(train.fig)
dev.off()

chunk.id <- basename(chunk.dir)
test.base <- file.path(chunks.dir, "figure-test-errors", chunk.id)
test.RData <- paste0(test.base, ".RData")
  
tobjs <- load(test.RData)

test.fig <- 
  ggplot()+
    ggtitle(animint.dir)+
    scale_y_continuous("aligned read coverage",
                       breaks=function(limits){
                         floor(limits[2])
                       })+
    scale_x_continuous(paste("position on", chrom,
                             "(kilo bases = kb)"))+
    coord_cartesian(xlim=lim.vec)+
    geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                      fill=annotation),
                  alpha=0.5,
                  color="grey",
                  data=regions)+
           scale_linetype_manual("error type",
                                 values=c(correct=0,
                                   "false negative"=3,
                                   "false positive"=1))+
    scale_fill_manual(values=ann.colors)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., labeller=function(var, val){
      sub("McGill0", "", sub(" ", "\n", val))
    }, scales="free")+
    geom_line(aes(base/1e3, count),
              data=some.counts,
              color="grey50")+
    geom_tallrect(aes(xmin=chromStart/1e3,
                      xmax=chromEnd/1e3,
                      linetype=status),
                  data=error.regions,
                  fill=NA,
                  size=1.5,
                  color="black")+
    geom_segment(aes(chromStart/1e3, 0,
                     xend=chromEnd/1e3, yend=0),
                 data=chunk.peaks, size=4, color="deepskyblue")+
    geom_point(aes(chromStart/1e3, 0),
               data=chunk.peaks, size=4, color="deepskyblue")

png.name <- paste0(test.base, ".png")
png(png.name, width=1000, h=(facet.rows+1)*30, units="px")
print(test.fig)
dev.off()

