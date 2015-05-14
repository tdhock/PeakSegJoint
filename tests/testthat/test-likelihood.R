library(testthat)
context("Poissson likelihood")

library(PeakSegJoint)

data(H3K36me3.TDH.other.chunk1)

some.counts <-
  subset(H3K36me3.TDH.other.chunk1$counts,
         43379893 < chromEnd)
profile.list <- ProfileList(some.counts)

test_that("C loss/segment means are correct", {
  param.values <- c(2, 3, 5)
  cDPA.seconds.list <- list()
  heuristic.seconds.list <- list()
  loss.list <- list()
  for(sample.i in seq_along(profile.list)){
    sample.id <- names(profile.list)[[sample.i]]
    cat(sprintf("%4d / %4d %s\n", sample.i, length(profile.list), sample.id))
    sample.counts <- profile.list[[sample.id]]
    one.profile.list <- profile.list[sample.id]
    m.args <- list()
    result.list <- list()
    for(param in param.values){
      param.name <- paste(param)
      param.int <- as.integer(param)
      fit <- PeakSegJointHeuristic(one.profile.list, param.int)
      model <- fit$models[[2]]
      expect_true(is.finite(model$loss))
      peak <- model$peak_start_end
      result.list[[param.name]] <- 
        data.frame(param.int,
                   chromStart=peak[1],
                   chromEnd=peak[2],
                   fit.loss=model$loss,
                   fit.seg1=model$seg1_mean_vec,
                   fit.seg2=model$seg2_mean_vec,
                   fit.seg3=model$seg3_mean_vec)
    }
    results <- do.call(rbind, result.list)

    count.vec <- with(sample.counts, rep(count, chromEnd-chromStart))
    sample.bases <- sum(with(sample.counts, chromEnd-chromStart))
    stopifnot(all.equal(sample.bases, length(count.vec)))
    seg1.chromStart <- sample.counts$chromStart[1]
    seg3.chromEnd <- sample.counts[nrow(sample.counts), "chromEnd"]
    seg1.chromEnd <- results$chromStart
    seg2.chromEnd <- results$chromEnd
    seg1.first <- seg1.chromStart - seg1.chromStart + 1
    seg1.last <- seg1.chromEnd - seg1.chromStart
    seg2.first <- seg1.last + 1
    seg2.last <- seg2.chromEnd - seg1.chromStart
    seg3.first <- seg2.last + 1
    seg3.last <- length(count.vec)
    data.sizes <- rowSums(cbind(seg3.last-seg2.last,
                                seg2.last-seg1.last,
                                seg1.last))
    stopifnot(data.sizes == length(count.vec))
    results$expected.loss <- results$expected.seg1 <-
      results$expected.seg2 <- results$expected.seg3 <- NA
    for(model.i in 1:nrow(results)){
      model <- results[model.i, ]
      seg.indices <-
        data.frame(first=c(seg1.first, seg2.first[model.i], seg3.first[model.i]),
                   last=c(seg1.last[model.i], seg2.last[model.i], seg3.last))
      seg.loss.list <- list()
      for(seg.i in 1:nrow(seg.indices)){
        col.name <- paste0("expected.seg", seg.i)
        s <- seg.indices[seg.i, ]
        seg.data <- count.vec[s$first:s$last]
        seg.mean <- results[[col.name]][model.i] <- mean(seg.data)
        seg.loss.list[[seg.i]] <- PoissonLoss(seg.data, seg.mean)
      }
      results$expected.loss[model.i] <- sum(unlist(seg.loss.list))
    }
    with(results, {
      expect_equal(fit.loss, expected.loss)
      expect_equal(fit.seg1, expected.seg1)
      expect_equal(fit.seg2, expected.seg2)
      expect_equal(fit.seg3, expected.seg3)
    })
  }
})

test_that("C flat models for 2 bin.factors agree", {
  fit2 <- PeakSegJointHeuristicStep1(profile.list, 2)
  fit3 <- PeakSegJointHeuristicStep1(profile.list, 3)
  expect_equal(fit2$sample_mean_vec, fit3$sample_mean_vec)
  expect_equal(fit2$flat_loss_vec, fit3$flat_loss_vec)

  for(profile.i in seq_along(profile.list)){
    one <- profile.list[profile.i]
    fit2 <- PeakSegJointHeuristicStep1(one, 2)
    fit3 <- PeakSegJointHeuristicStep1(one, 3)
    expect_equal(fit2$sample_mean_vec, fit3$sample_mean_vec)
    expect_equal(fit2$flat_loss_vec, fit3$flat_loss_vec)
    pro <- profile.list[[profile.i]]
    bases <- as.numeric(with(pro, chromEnd-chromStart))

    mean.value <- sum(as.numeric(pro$count)*bases)/sum(bases)
    expect_equal(fit2$sample_mean_vec, mean.value)
    expect_equal(fit3$sample_mean_vec, mean.value)

    loss.value <- PoissonLoss(pro$count, mean.value, bases)
    expect_equal(fit2$flat_loss_vec, loss.value)
    expect_equal(fit3$flat_loss_vec, loss.value)
  }
})
