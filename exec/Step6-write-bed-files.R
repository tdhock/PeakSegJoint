library(PeakSegJoint)
library(data.table)
library(animint)

argv <- system.file("exampleData", "PeakSegJoint-predictions",
                    package="PeakSegJoint")

argv <- "~/genomepipelines/H3K4me3_TDH_immune/PeakSegJoint-predictions"
argv <- "~/projects/blueprint/results/H3K27ac/PeakSegJoint-predictions"
argv <- "~/projects/input-test-data-master/PeakSegJoint-predictions/"

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

if(length(argv) != 1){
  stop("usage: Step6.R path/to/PeakSegJoint-predictions")
}

pred.dir <- normalizePath(argv[1], mustWork=TRUE)
data.dir <- dirname(pred.dir)

positive.regions.RData <- file.path(data.dir, "positive.regions.RData")
load(positive.regions.RData)

RData.file.vec <- Sys.glob(file.path(pred.dir, "*.RData"))

peak.mat.list <- list()
pred.peak.list <- list()
for(RData.file in RData.file.vec){
  objs <- load(RData.file)
  pred.peak.list[[RData.file]] <- data.table(pred.peaks)
  peak.mat.list[[RData.file]] <- peak.mat
}
all.peaks.mat <- do.call(cbind, peak.mat.list)
all.peaks.dt <- do.call(rbind, pred.peak.list)
all.peaks.dt[, sample.path := file.path(sample.group, sample.id)]
all.peaks.df <- data.frame(all.peaks.dt)

## Join positive region labels with peak calls.
Input.counts <- all.peaks.dt[, list(
  Input=sum(sample.group=="Input")
  ), by=.(peak.name, chrom, chromStart, chromEnd)]
if(is.data.frame(positive.regions)){
  setkey(Input.counts, chrom, chromStart, chromEnd)
  positive.dt <- data.table(positive.regions)
  setkey(positive.dt, chrom, regionStart, regionEnd)
  over.dt <- foverlaps(Input.counts, positive.dt, nomatch=0L)
  thresh.dt <- over.dt[, list(
    specific=sum(!input.has.peak),
    nonspecific=sum(input.has.peak)
    ), by=Input]
  ## The decision function calls a peak "specific" if at most T input
  ## samples have peaks, where T is the learned threshold: in the error
  ## rate table below, [min.thresh, max.thresh), e.g. [4, 0) means that
  ## thresholds 4, 3, 2, 1 all have the same error. 
  specific.error.dt <- thresh.dt[, {
    min.thresh <- c(-Inf, Input)
    max.thresh <- c(Input, Inf)
    threshold <- floor((min.thresh+max.thresh)/2)
    threshold[1] <- min(Input)
    threshold[length(threshold)] <- max(Input)
    data.table(
      min.thresh,
      max.thresh,
      threshold,
      FN=c(0, cumsum(rev(specific))),
      FP=c(rev(cumsum(nonspecific)), 0))
  }]
  specific.error.dt[, errors := FP + FN]
  max.input.samples <- specific.error.dt[which.min(errors), threshold]
  thresh.dt <- data.table(max.input.samples, thresh.type="specific")
  specific.error <- data.frame(specific.error.dt)
  print(specific.error)
  cat('A peak is called "specific" and included in output files if at most',
      max.input.samples,
      "Input samples have peaks. Counts:\n")
  Input.counts[, call := ifelse(Input <= max.input.samples,
                        "specific", "nonspecific")]
  print(table(Input.counts$call))
  specific.peak.names <- Input.counts[call=="specific", peak.name]
}else{
  cat("All peaks included in output files,",
      "since there are no labeled peaks in Input samples.\n")
  specific.peak.names <- colnames(all.peaks.mat)
  specific.error <- NULL
}

writeBoth <- function(peaks.df, bed.path, header.line){
  sorted.peaks <- data.table(peaks.df)[order(chrom, chromStart), ]
  bigBed.file <- sub("[.]bed$", ".bigBed", bed.path)
  write.table(sorted.peaks, bed.path, 
              quote=FALSE, row.names=FALSE, col.names=FALSE)
  s <- ifelse(nrow(sorted.peaks)==1, "", "s")
  message("wrote ", nrow(sorted.peaks), " peak", s, " to ", bed.path)
  bigBedCmd <- paste("bedToBigBed", bed.path, chrom.file, bigBed.file)
  status <- system(bigBedCmd)
  if(status != 0){
    stop("error in bedToBigBed")
  }
  ## Also write a bed.gz file with a header for uploading as a custom
  ## track.
  bed.gz <- paste0(bed.path, ".gz")
  con <- gzfile(bed.gz, "w")
  writeLines(header.line, con)
  write.table(sorted.peaks, con, 
              quote=FALSE, row.names=FALSE, col.names=FALSE)
  close(con)
}
bigwig.file.vec <- Sys.glob(file.path(data.dir, "*", "*.bigwig"))

chrom.ranges <- bigWigInfo(bigwig.file.vec[1])
chrom.vals <- unique(c(
  paste0("chr", c(1:22, "X", "Y", "M")),
  chrom.ranges$chrom))
chrom.ranges[, chrom := factor(chrom, chrom.vals)]
just.ends <- data.frame(chrom.ranges)[, c("chrom", "chromEnd")]
chrom.file <- tempfile()
write.table(just.ends, chrom.file, quote=FALSE,
            row.names=FALSE, col.names=FALSE)

## Viz for each cell type.
sample.type <- factor(sub("[0-9]*/.*", "", rownames(all.peaks.mat)))
sample.counts <- table(sample.type)
sample.counts.dt <- data.table(sample.counts, thresh.type="max samples")
counts.not.Input <- sample.counts.dt[sample.type != "Input", ]
counts.not.Input[, nonInputType := sample.type]
counts.Input <- sample.counts.dt[sample.type == "Input", ]
type.mat <- matrix(sample.type, nrow(all.peaks.mat), ncol(all.peaks.mat),
                   dimnames=dimnames(all.peaks.mat))
type.mat[all.peaks.mat==0] <- NA
type.count.mat <- apply(type.mat, 2, function(x){
  table(factor(x, levels(sample.type)))
})
u.peaks.df <-
  unique(all.peaks.df[, c("chrom", "chromStart", "chromEnd", "peak.name")])
rownames(u.peaks.df) <- u.peaks.df$peak.name
u.peaks.df <- u.peaks.df[colnames(all.peaks.mat), ]
stopifnot(u.peaks.df$peak.name == colnames(all.peaks.mat))
stopifnot(!is.na(u.peaks.df))
count.dt.list <- list()
for(type in rownames(type.count.mat)){
  count.dt.list[[type]] <- 
    data.table(type,
               u.peaks.df,
               samples.up=type.count.mat[type, ])
}
count.dt <- do.call(rbind, count.dt.list)[0 < samples.up, ]
type.levels <- unique(c("Input", sample.counts.dt$sample.type))
count.dt[, type.fac := factor(type, type.levels)]
setkey(count.dt, peak.name)

type.vec <- rownames(type.count.mat)[rownames(type.count.mat)!="Input"]
scatter.dt.list <- list()
chromCounts.list <- list()
for(nonInputType in type.vec){
  type.dt <- 
    data.table(up=type.count.mat[nonInputType, ],
               peak.name=colnames(type.count.mat))
  by.vec <- if("Input" %in% rownames(type.count.mat)){
    type.dt$Input <- type.count.mat["Input", ]
    type.dt[, dotID := sprintf("%d %s samples, %d Input samples",
                      up, nonInputType, Input)]
    c("dotID", "up", "Input")
  }else{
    type.dt[, dotID := sprintf("%d %s samples", up, nonInputType)]
    c("dotID", "up")
  }
  some.up <- type.dt[0 < up, ]
  setkey(some.up, peak.name)
  some.up.counts <- some.up[, list(count=.N), by=by.vec]
  scatter.dt.list[[nonInputType]] <- data.table(nonInputType, some.up.counts)
  chromCounts.list[[nonInputType]] <-
    data.table(nonInputType, count.dt[some.up])
}
chromCounts <- do.call(rbind, chromCounts.list)
scatter.dt <- do.call(rbind, scatter.dt.list)

scatter.dt[, totals := paste(count, "regions with a common peak in", dotID)]
scatter.dt[, chrom := factor("chrY", chrom.vals)]
chrom.ranges[, ymin := factor(type.levels[1], type.levels) ]
chrom.ranges[, ymax := factor(type.levels[length(type.levels)], type.levels) ]
chromCounts[, bases := chromEnd-chromStart]
chromCounts[, middle := as.integer((chromStart+chromEnd)/2)]
chromCounts[, zoom.bases := bases * 10]
chromCounts[, zoomStart := as.integer(middle-zoom.bases)]
chromCounts[, zoomEnd := as.integer(middle+zoom.bases)]
setkey(chrom.ranges, chrom)
chrom.size.vec <- chrom.ranges[paste(chromCounts$chrom), chromEnd]
chromCounts[, relative.middle := middle / chrom.size.vec]
chromCounts[, chrom := factor(chrom, chrom.vals)]
countsByChrom <- chromCounts[, list(peaks=length(unique(peak.name))),
                             by=.(dotID, chrom)]

setkey(countsByChrom, chrom)
bg.rect.list <- list()
for(chrom in chrom.vals){
  text.for.chrom <- countsByChrom[chrom]
  rect.for.chrom <- data.table(chrom, scatter.dt)
  rect.for.chrom[, size := ifelse(dotID %in% text.for.chrom$dotID, 0.55, 0.4)]
  bg.rect.list[[chrom]] <- rect.for.chrom
}
bg.rect <- do.call(rbind, bg.rect.list)

## PredictedPeaks <- list(
##   chromCounts=data.frame(chromCounts),
##   countsByChrom=data.frame(countsByChrom),
##   chrom.ranges=data.frame(chrom.ranges),
##   scatter.text=data.frame(scatter.dt),
##   counts.Input=data.frame(counts.Input),
##   counts.not.Input=data.frame(counts.not.Input),
##   bg.rect=data.frame(bg.rect))
## save(PredictedPeaks, file="~/R/animint-examples/data/PredictedPeaks.RData")

hgTracks <- "http://genome.ucsc.edu/cgi-bin/hgTracks"
viz <- list(
  title=pred.dir,
  oneChrom=ggplot()+
    ggtitle("PeakSegJoint detections on selected chromosome")+
    theme_bw()+
    coord_cartesian(xlim=c(0, 1))+
    theme_animint(width=1500,
                  height=50 + 20 * length(unique(chromCounts$type)))+
    theme(axis.line.x=element_blank(), axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), axis.title.x=element_blank())+
    scale_y_discrete("cell type", drop=FALSE)+
    geom_text(aes(relative.middle, type.fac, label=samples.up,
                  href=sprintf("%s?db=hg19&position=%s:%d-%d",
                    hgTracks, chrom, zoomStart, zoomEnd),
                  showSelected2=chrom,
                  showSelected=dotID),
              size=11,
              data=chromCounts),
  ## TODO: rather than displaying the entire chromosome on the top
  ## plot (which often squishes together a bunch of nearby peak
  ## counts), we can tile the genome with boxes as in
  ## https://github.com/tdhock/animint-examples/blob/master/examples/scaffolds.R
  chroms=ggplot()+
    theme_bw()+
    theme_animint(width=1500, height=330)+
    scale_y_discrete("chromosome", drop=FALSE)+ 
    scale_x_continuous("position on chromosome (mega bases)")+
    geom_text(aes(0, chrom, label=paste0(peaks, "_"),
                  clickSelects=chrom,
                  showSelected=dotID),
              hjust=1,
              size=11,
              data=countsByChrom)+
    geom_segment(aes(chromStart/1e6, chrom,
                     clickSelects=chrom,
                     xend=chromEnd/1e6, yend=chrom),
              size=9,
              data=chrom.ranges)+
    geom_point(aes(chromEnd/1e6, chrom,
                   clickSelects=chrom),
              size=5,
              data=chrom.ranges)+
    geom_text(aes(max(chrom.ranges$chromEnd)/2e6, chrom,
                  showSelected=dotID,
                  label=totals),
             data=scatter.dt))

if("Input" %in% rownames(type.count.mat)){
  p <- ggplot()+
    geom_vline(aes(xintercept=N, color=thresh.type),
               data=counts.not.Input)+
    geom_hline(aes(yintercept=N, color=thresh.type),
               show_guide=TRUE,
               data=counts.Input)
  if(is.data.frame(positive.regions)){
    p <- p+
      geom_hline(aes(yintercept=max.input.samples+0.5, color=thresh.type),
                 show_guide=TRUE,
                 data=thresh.dt)
  }
  p <- p+
    scale_color_manual("threshold", values=c(
                                      "max samples"="grey",
                                      specific="grey30"))+
    scale_x_continuous("number of samples with a peak")+
    facet_grid(nonInputType ~ .)+
    theme_bw()+
    scale_fill_gradient(low="grey", high="red")+
    theme_animint(width=1500,
                  height=50 + 100 * length(unique(chromCounts$nonInputType)))+
    theme(panel.margin=grid::unit(0, "cm"))+
    geom_rect(aes(xmin=up-size, xmax=up+size,
                  ymin=Input-size, ymax=Input+size,
                  tooltip=totals,
                  clickSelects=dotID,
                  showSelected=chrom,
                  fill=log10(count)),
              color="transparent",
              data=bg.rect)
  viz$scatter <- p
}else{
  viz$bars <- ggplot()+
    geom_vline(aes(xintercept=N, color=thresh.type),
               show_guide=TRUE,
               data=counts.not.Input)+
    scale_color_manual("threshold", values=c(
                                      "max samples"="grey",
                                      specific="grey30"))+
    scale_x_continuous("number of samples with a peak")+
    scale_y_continuous("count")+
    facet_grid(nonInputType ~ .)+
    theme_bw()+
    theme_animint(width=1500,
                  height=50 + 100 * length(unique(chromCounts$nonInputType)))+
    theme(panel.margin=grid::unit(0, "cm"))+
    geom_bar(aes(up, count,
                 tooltip=totals,
                 clickSelects=dotID),
             stat="identity",
             data=scatter.dt)
}  

viz.dir <- paste0(pred.dir, "-viz")
animint2dir(viz, viz.dir)

## Save a summary bed file for easy viewing of the entire model as a
## custom track on UCSC.
specific.peaks.df <- subset(all.peaks.df, peak.name %in% specific.peak.names)
u.cols <- c("chrom", "chromStart", "chromEnd")
u.peaks <-
  unique(specific.peaks.df[, u.cols])
rownames(u.peaks) <- with(u.peaks, {
  paste0(chrom, ":", chromStart, "-", chromEnd)
})
specific.peaks.df$group <- sub("[0-9]*$", "", specific.peaks.df$sample.group)
count.mat <- with(specific.peaks.df, table(group, peak.name))
getCount <- function(x){
  count.str <- paste0(x, ":", rownames(count.mat))
  non.zero <- count.str[x != 0]
  paste(non.zero, collapse="/")
}
str.vec <- apply(count.mat, 2, getCount)
u.peaks$name <- str.vec[rownames(u.peaks)]
count.tab <- colSums(all.peaks.mat)
u.peaks$score <- count.tab[rownames(u.peaks)]
summary.bed <- file.path(data.dir, "PeakSegJoint.summary.bed")
header <- 
  paste("track",
        "visibility=pack",
        "name=PeakSegJointGroups",
        'description="PeakSegJoint sample counts by group"')
writeBoth(u.peaks, summary.bed, header)

predictions.RData <- file.path(data.dir, "PeakSegJoint.predictions.RData")
predictions.csv <- file.path(data.dir, "PeakSegJoint.predictions.csv")
peak.Input.counts <- data.frame(Input.counts)
save(all.peaks.mat, all.peaks.df,
     specific.peak.names,
     peak.Input.counts, specific.error,
     file=predictions.RData)
specific.peaks.mat <- all.peaks.mat[, specific.peak.names]
write.csv(specific.peaks.mat, predictions.csv, quote=TRUE,
          row.names=TRUE)

peaks.by.path <- split(specific.peaks.df, specific.peaks.df$sample.path)
for(sample.path in names(peaks.by.path)){
  sample.peaks <- peaks.by.path[[sample.path]]
  ord <- with(sample.peaks, order(chrom, chromStart))
  sorted.peaks <- sample.peaks[ord, ]
  bed.peaks <- sorted.peaks[, c("chrom", "chromStart", "chromEnd")]
  rownames(bed.peaks) <- NULL
  bed.file <- file.path(data.dir, paste0(sample.path, ".bed"))
  desc.line <-
    sprintf('description="PeakSegJoint calls for sample %s"',
            sample.path)
  header <- 
    paste("track",
          "visibility=pack",
          paste0("name=PeakSegJoint:", sample.path),
          desc.line)
  writeBoth(bed.peaks, bed.file, header)
}
