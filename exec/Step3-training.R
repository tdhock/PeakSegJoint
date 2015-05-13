require(data.table)
require(ggplot2)
require(xtable)
require(PeakSegJoint)

argv <-
  system.file(file.path("exampleData",
                        "PeakSegJoint-chunks"),
              package="PeakSegDP")

argv <- "PeakSegJoint-chunks/H3K36me3_TDH_immune"

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

if(length(argv) != 1){
  stop("usage: Step2.R path/to/PeakSegJoint-chunks")
}

chunks.dir <- normalizePath(argv[1])

problems.RData.vec <- Sys.glob(file.path(chunks.dir, "*", "problems.RData"))

err.mat.list <- list()
chunk.best.list <- list()
bpp.list <- list()
for(chunk.i in seq_along(problems.RData.vec)){
  problems.RData <- problems.RData.vec[[chunk.i]]
  chunk.dir <- dirname(problems.RData)
  chunk.id <- basename(chunk.dir)
  index.html <-
    file.path("..", chunk.id, "figure-train-errors", "index.html")
  objs <- load(problems.RData)
  err.vec <- step2.error$errors
  min.df <- subset(step2.error, errors==min(errors))
  min.errors <- min(step2.error$errors)
  bpp.list[[chunk.id]] <- min.df$bases.per.problem
  href <- sprintf('<a href="%s">%s</a>', index.html, chunk.id)
  chunk.best.list[[chunk.id]] <-
    data.frame(chunk.id=href,
               min.errors,
               regions=step2.error$regions[1])
  names(err.vec) <- rownames(step2.error)
  err.mat.list[[chunk.i]] <- err.vec
}
chunk.best <- do.call(rbind, chunk.best.list)
err.mat <- do.call(rbind, err.mat.list)
err.vec <- colSums(err.mat)
err.df <- data.frame(bases.per.problem=as.integer(names(err.vec)),
                     errors=err.vec)
chosen.i <- which.min(err.vec)
res.str <- names(err.vec)[chosen.i]
chosen.df <- err.df[chosen.i, ]

chunk.best$selected.errors <- err.mat[, res.str]
chunk.best$bases.per.problem <- NA
for(chunk.i in seq_along(bpp.list)){
  bpp <- bpp.list[[chunk.i]]
  bpp[bpp == res.str] <-
    paste0("<b>", bpp[bpp==res.str], "</b>")
  chunk.best$bases.per.problem[chunk.i] <- paste0(bpp, collapse="<br />")
}

chunk.ordered <-
  chunk.best[order(chunk.best$selected.errors, decreasing=TRUE),]

resCurve <- 
ggplot()+
  geom_line(aes(bases.per.problem, errors),
            data=err.df)+
  scale_x_log10()+
  ggtitle(paste(res.str, "bases/problem"))+
  geom_point(aes(bases.per.problem, errors),
             data=chosen.df,
             pch=1)

png.name <-
  file.path(chunks.dir, "figure-train-errors", "figure-train-errors.png")
png.dir <- dirname(png.name)
dir.create(png.dir, showWarnings=FALSE, recursive=TRUE)
png(png.name, width=400, h=200, units="px")
print(resCurve)
dev.off()

res.xt <- xtable(err.df)
res.html <-
  print(res.xt, type="html",
        include.rownames=FALSE)

xt <- xtable(chunk.ordered)
html.table <-
  print(xt, type="html",
        include.rownames=FALSE,
        sanitize.text.function=identity)
html.out <-
  paste(res.html,
        '<br /> <img src="figure-train-errors.png" /> <br />',
        html.table)

out.file <- file.path(chunks.dir, "figure-train-errors", "index.html")

cat(html.out, file=out.file)

problems.by.chunk <- list()
for(chunk.i in seq_along(problems.RData.vec)){
  problems.RData <- problems.RData.vec[[chunk.i]]
  objs <- load(problems.RData)
  chunk.problems <- step2.by.res[[res.str]]
  for(problem.name in names(chunk.problems)){
    info <- chunk.problems[[problem.name]]
    target <- info$error$target
    n.finite <- sum(is.finite(target))
    if(n.finite > 0){
      problems.by.chunk[[paste(chunk.i)]][[problem.name]] <-
        list(features=info$features,
             target=target)
    }
  }
}

print("TODO IntervalRegressionProblems")
