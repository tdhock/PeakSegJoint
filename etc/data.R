load("~/projects/chip-seq-paper/chunks/H3K36me3_TDH_other/1/counts.RData")
load("~/projects/chip-seq-paper/chunks/H3K36me3_TDH_other/1/regions.RData")

some.counts <- subset(counts, 43100000 < chromEnd & chromStart < 43200000)
some.regions <- subset(regions, chromStart < 43200000)

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

ggplot()+
  scale_fill_manual(values=ann.colors)+
  geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                    fill=annotation),
                data=some.regions,
                alpha=0.5)+
  geom_step(aes(chromStart/1e3, coverage), data=some.counts)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(sample.id ~ ., scales="free")

counts$count <- counts$coverage
meta.cols <- c("cell.type", "sample.id", "chromStart", "chromEnd")
H3K36me3.TDH.other.chunk1 <-
  list(counts=counts[,c(meta.cols, "count")],
       regions=regions[, c(meta.cols, "annotation")])

save(H3K36me3.TDH.other.chunk1, file="../data/H3K36me3.TDH.other.chunk1.RData")
