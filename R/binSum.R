binSum <- structure(function
### Compute sum of compressed coverage profile in bins, using fast C
### code.
(compressed,
### data.frame with integer columns chromStart, chromEnd, count.
 bin.chromStart=0L,
### Base before first bin.
 bin.size=1L,
### Bin size.
 n.bins=2000L,
### Number of bins.
 empty.as.zero=FALSE
### Sometimes the last few bins do not have any overlapping data in
### compressed. If TRUE, set these counts to 0. If FALSE, ignore these
### bins (returning a data.frame with fewer than n.bins rows).
 ){
  stopifnot(is.integer(compressed$chromStart))
  stopifnot(is.integer(compressed$chromEnd))
  stopifnot(is.integer(compressed$count))
  stopifnot(is.integer(bin.chromStart))
  stopifnot(length(bin.chromStart) == 1)
  stopifnot(is.integer(bin.size))
  stopifnot(length(bin.size) == 1)
  stopifnot(is.integer(n.bins))
  stopifnot(length(n.bins) == 1)
  bin.total <- integer(n.bins)
  result <- 
  .C("binSum_interface",
     profile.chromStart=as.integer(compressed$chromStart),
     profile.chromEnd=as.integer(compressed$chromEnd),
     profile.coverage=as.integer(compressed$count),
     n.profiles=as.integer(nrow(compressed)),
     bin.total=as.integer(bin.total),
     bin.size=as.integer(bin.size),
     n.bins=as.integer(n.bins),
     bin.chromStart=as.integer(bin.chromStart),
     package="PeakSegDP")
  total <- result$bin.total
  if(empty.as.zero){
    total[total == -1] <- 0L
  }
  chromStart <- seq(bin.chromStart, by=bin.size, l=n.bins)
  chromEnd <- chromStart + bin.size
  data.frame(chromStart, chromEnd, count=total,
             mean=total/bin.size)[total >= 0, ]
### data.frame with n.bins rows and columns chromStart, chromEnd,
### count, mean.
}, ex=function(){
  ## bins of size 3bp.
  ## -1-   -3-   -5-
  ##    -2-   -4-
  ## 123456789012345 base index.
  ## --2---
  ##       --1-
  ##           --0-------
  ## Coverage profile.
  profile <- data.frame(chromStart=as.integer(c(0, 6, 10)),
                        chromEnd=as.integer(c(6, 10, 10000)),
                        count=as.integer(c(2, 1, 0)))
  library(PeakSegJoint)
  bins <- binSum(profile,
                 bin.chromStart=0L,
                 bin.size=3L,
                 n.bins=2000L)
  library(ggplot2)
  bases <- data.frame(position=1:15, base="N")
  ggplot()+
    ylab("")+
    geom_text(aes(position, 0, label=base),
              data=bases)+
    geom_step(aes(chromStart+0.5, count, color=what),
              data=data.frame(profile, what="profile"),
              size=2)+
    geom_step(aes(chromStart+0.5, count, color=what),
              data=data.frame(bins, what="bin total"))+
    geom_step(aes(chromStart+0.5, mean, color=what),
              data=data.frame(bins, what="bin mean"))+
    coord_cartesian(xlim=c(0, max(bases$position)))
})
