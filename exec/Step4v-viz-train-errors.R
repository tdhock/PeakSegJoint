library(animint)

## Filter counts by limits (area around regions).
limits <- regions[, .(chromStart=min(chromStart),
                      chromEnd=max(chromEnd))]
limits[, bases := chromEnd - chromStart ]
limits[, expand := as.integer(bases/10) ]
limits[, min := chromStart - expand]
limits[, max := chromEnd + expand]
lim.vec <- with(limits, c(min, max))/1e3
some.counts <- counts[limits$min < chromEnd &
                        chromStart < limits$max,]

viz <-
  list(coverage=ggplot()+
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
         theme_animint(width=1000, height=(length(counts.by.sample)+2)*100)+
         theme(panel.margin=grid::unit(0, "cm"))+
         facet_grid(sample.id ~ ., labeller=function(var, val){
           sub("McGill0", "", sub(" ", "\n", val))
         }, scales="free")+
         geom_step(aes(chromStart/1e3, count),
                   data=some.counts,
                   color="grey50"))


dvec <- diff(log(res.error$bases.per.problem))
dval <- exp(mean(dvec))
dval2 <- (dval-1)/2 + 1
viz$resError <- ggplot()+
  geom_tallrect(aes(xmin=bases.per.problem/dval2,
                    xmax=bases.per.problem*dval2,
                    clickSelects=bases.per.problem),
             alpha=0.5,
             data=res.error)+
  scale_x_log10()+
  geom_line(aes(bases.per.problem, errors),
            data=res.error)
viz$title <- chunk.dir

animint.dir <- file.path(chunk.dir, "figure-train-errors")
animint2dir(viz, animint.dir, open.browser=FALSE)

