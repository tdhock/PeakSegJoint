library(animint)
library(PeakSegJoint)

argv <- c(
  "~/projects/blueprint/results/H3K27ac/PeakSegJoint.predictions.RData",
  "chr5:156569171-156570800")

argv <- c(
  "~/input-test-data-processed/PeakSegJoint.predictions.RData",
  "chr10:33183540-33225283")

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

if(length(argv) != 2){
  stop("usage: PlotPeak.R PeakSegJoint.predictions.RData chrX:1234-5678")
}

facet.labels <- function(var, val){
  sub("_.*", "", val)
}

predictions.RData <- normalizePath(argv[1], mustWork=TRUE)
peak.to.plot <- argv[2]
data.dir <- dirname(predictions.RData)

load(predictions.RData)

stopifnot(peak.to.plot %in% colnames(all.peaks.mat))
center.peaks <- subset(all.peaks.df, peak.name==peak.to.plot)

one.peak <- with(center.peaks[1,], {
  data.table(chrom, peakStart=chromStart, peakEnd=chromEnd)
})
zoom.bases <- one.peak[, peakEnd - peakStart] * 1
one.peak[, zoomStart := peakStart - zoom.bases]
one.peak[, zoomEnd := peakEnd + zoom.bases]

## Should also show other peaks in the zoomed-out region.
all.peaks.dt <- data.table(all.peaks.df)
setkey(all.peaks.dt, chrom, chromStart, chromEnd)
setkey(one.peak, chrom, zoomStart, zoomEnd)
peaks.to.plot <- foverlaps(all.peaks.dt, one.peak, nomatch=0L)

bigwig.glob <- file.path(data.dir, "*", "*.bigwig")
bigwig.file.vec <- Sys.glob(bigwig.glob)
counts.by.sample <- list()
for(bigwig.i in seq_along(bigwig.file.vec)){
  bigwig.file <- bigwig.file.vec[[bigwig.i]]
  cat(sprintf("%4d / %4d bigwigs %s\n", bigwig.i, length(bigwig.file.vec),
              bigwig.file))
  sample.counts <- one.peak[, {
    readBigWig(bigwig.file, chrom, zoomStart, zoomEnd)
  }]
  if(nrow(sample.counts)){
    sample.id <- sub("[.]bigwig$", "", basename(bigwig.file))
    sample.group <- basename(dirname(bigwig.file))
    counts.by.sample[[paste(sample.id, sample.group)]] <-
      data.table(sample.id, sample.group, sample.counts)
  }
}
some.counts <- do.call(rbind, counts.by.sample)
if(is.null(some.counts)){
  stop("no coverage data to plot in ", bigwig.glob)
}

n.profiles <- length(counts.by.sample)
height.pixels <- (n.profiles+1)*30

gg <- 
  ggplot()+
    ggtitle(paste("Peak", peak.to.plot,
                  "predicted in", nrow(center.peaks),
                  "samples"))+
    scale_y_continuous("aligned read coverage",
                       breaks=function(limits){
                         floor(limits[2])
                       })+
    scale_x_continuous(paste("position on", one.peak$chrom,
                             "(kilo bases = kb)"))+
    coord_cartesian(xlim=one.peak[, c(zoomStart, zoomEnd)/1e3])+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.group + sample.id ~ ., labeller=facet.labels,
               scales="free")+
    geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                  ymin=0, ymax=count),
              color="grey50",
              data=some.counts)+
    geom_segment(aes(chromStart/1e3, 0,
                     xend=chromEnd/1e3, yend=0),
                 data=peaks.to.plot,
                 size=4,
                 color="deepskyblue")

peaksPlot <- list(
  peaks=peaks.to.plot,
  coverage=some.counts)

out.dir <- file.path(data.dir, "PeakSegJoint-predictions-plots", peak.to.plot)
dir.create(out.dir, showWarnings=FALSE, recursive=TRUE)

peaksPlot.RData <- file.path(out.dir, "peaksPlot.RData")
save(peaksPlot, file=peaksPlot.RData)

png.name <- file.path(out.dir, "peaksPlot.png")
png(png.name, width=1000, h=height.pixels, units="px")
print(gg)
dev.off()

