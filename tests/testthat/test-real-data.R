context("real data")

library(PeakSegJoint)

data(H3K36me3.AM.immune.chunk21)

test_that("21 peak loss < 20 peak loss", {
  fit <- PeakSegJointHeuristic(H3K36me3.AM.immune.chunk21)
  converted <- ConvertModelList(fit)
  with(converted$loss, {
    expect_that(loss[22], is_less_than(loss[21]))
  })

  segs <- subset(converted$segments, peaks %in% c(20, 21))
  segs$peaks.str <- paste(segs$peaks)

  ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., labeller=function(var, val){
      sub("McGill0", "", sub(" ", "\n", val))
    }, scales="free")+
    geom_step(aes(chromStart/1e3, count),
              data=H3K36me3.AM.immune.chunk21,
              color="grey50")+
    geom_segment(aes(chromStart/1e3, mean,
                     xend=chromEnd/1e3, yend=mean,
                     size=peaks.str,
                     color=peaks.str),
                 data=segs)+
    scale_size_manual(values=c("21"=1, "20"=2))
  
})

data(H3K4me3.TDH.other.chunk8)

test_that("8 peaks are feasible", {
  fit <- PeakSegJointHeuristic(H3K4me3.TDH.other.chunk8)
  converted <- ConvertModelList(fit)

  peaks <- unique(converted$peaks[, c("peaks", "chromStart", "chromEnd")])

  peaks7 <- subset(converted$peaks, peaks == 7)
  expect_false("McGill0267" %in% peaks7$sample.id)

  segs <- subset(converted$segments, peaks == 7)
  breaks <- subset(segs, chromStart > min(chromStart))
  
  ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., labeller=function(var, val){
      sub("McGill0", "", sub(" ", "\n", val))
    }, scales="free")+
    geom_step(aes(chromStart/1e3, count),
              data=H3K4me3.TDH.other.chunk8,
              color="grey50")+
    geom_vline(aes(xintercept=chromStart/1e3),
               linetype="dashed",
               color="green",
               data=breaks)+
    geom_segment(aes(chromStart/1e3, mean,
                     xend=chromEnd/1e3, yend=mean),
                 color="green",
                 size=1,
                 data=segs)
})
