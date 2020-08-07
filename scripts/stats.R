library(methylKit)

setwd('data')

list.files = list('data/sorted_index/EV1_sorted.bam', 'data/sorted_index/set1rep1_sorted.bam')

objs = processBismarkAln(
        list.files,
	sample.id = list('ev1', 'setrep1'),
	treatment = c(1,0),
	read.context = 'CpG')

pdf('methylationStats.pdf')

for (i in 1:length(list.files)){
    getMethylationStats(objs[[i]], plot=T, both.strands=F}
