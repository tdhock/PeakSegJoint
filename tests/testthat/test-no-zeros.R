library(testthat)
context("sparse data without counts equal to 0")

library(PeakSegJoint)

data(peak1.infeasible)
data.list <-
  list(with.zero=peak1.infeasible,
       without.zero=subset(peak1.infeasible, 0 < count))

test_that("binSum gives same answer for real data sets", {
  bin.chromStart <- min(peak1.infeasible$chromStart)
  bin.chromEnd <- max(peak1.infeasible$chromEnd)
  bases.per.bin <- 100L
  n.bins <- as.integer((bin.chromEnd-bin.chromStart)/bases.per.bin)
  data.sample.list <- lapply(data.list, function(df)split(df, df$sample.id))
  sample.id.vec <- names(data.sample.list[[1]])
  for(sample.id in sample.id.vec){
    bin.list <- list()
    sample.list <- list()
    for(data.name in names(data.sample.list)){
      one.sample <- data.sample.list[[data.name]][[sample.id]]
      bins <- binSum(one.sample, bin.chromStart, bases.per.bin, n.bins)
      bin.list[[data.name]] <- data.frame(bins, data.name)
      sample.list[[data.name]] <- data.frame(one.sample, data.name)
    }
    for(data.name in names(data.sample.list)){
      bin.list[[data.name]]$difference <-
        ifelse(bin.list[[data.name]]$count != bin.list$with.zero$count,
               "different", "same")
    }
    both.samples <- do.call(rbind, sample.list)
    both.bins <- do.call(rbind, bin.list)
    head(bin.list$without.zero, 20)
    ggplot()+
      ##coord_cartesian(xlim=c(118219186, 118219286)/1e3)+
      geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                    ymin=0, ymax=count),
                data=both.samples)+
      geom_step(aes(chromStart/1e3, count),
                data=both.samples)+
      geom_segment(aes(chromStart/1e3, mean,
                       color=difference,
                       xend=chromEnd/1e3, yend=mean),
                   data=both.bins)+
      theme_bw()+
      theme(panel.margin=grid::unit(0, "cm"))+
      facet_grid(data.name ~ .)
    with(bin.list, rbind(with.zero$count, without.zero$count))
    with(bin.list, expect_equal(with.zero$count, without.zero$count))
  }
})

test_that("PeakSegJointHeuristic loss same with or without zeros", {
  for(bp in 2:7){
    loss.list <- list()
    for(data.name in names(data.list)){
      L <- data.list[[data.name]]
      fit <- PeakSegJointHeuristic(L, bp)
      loss.list[[data.name]] <- sapply(fit$models, "[[", "loss")
    }
    with(loss.list, expect_equal(with.zero, without.zero))
  }
})


test_that("PeakSegJointHeuristicStep1 loss same with or without zeros", {
  for(bp in 2:7){
    loss.list <- list()
    for(data.name in names(data.list)){
      L <- data.list[[data.name]]
      fit <- PeakSegJointHeuristicStep1(L, bp)
      loss.list[[data.name]] <- sapply(fit$models, "[[", "loss")
    }
    with(loss.list, expect_equal(with.zero, without.zero))
  }
})

bp <- 4
peak.list <- list()
for(data.name in names(data.list)){
  L <- data.list[[data.name]]
  fit <- PeakSegJointHeuristic(L, bp)
  peak.list[[data.name]] <- sapply(fit$models, "[[", "peak_start_end")
}

test_that("PeakSegJointSeveral loss same with or without zeros", {
  loss.list <- list()
  for(data.name in names(data.list)){
    L <- data.list[[data.name]]
    fit <- PeakSegJointSeveral(L)
    loss.vec <- sapply(fit$models, "[[", "loss")
    expect_true(all(is.finite(loss.vec)))
    loss.list[[data.name]] <- loss.vec
  }
  with(loss.list, expect_equal(with.zero, without.zero))
})

