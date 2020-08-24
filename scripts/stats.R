if(!require('methylKit')){
    library('BiocManager')
    BiocManager::install('methylKit')
    library('methylKit')
}
if(!require('genomation')){
    BiocManager::install('genomation')
    library('genomation')
}

#setwd(paste0(getwd(), '/data/output'))

list.files = list('data/sorted_index/EV1_sorted.bam', 'data/sorted_index/set1rep1_sorted.bam', 'data/sorted_index/by4741_sorted.bam')

objs = processBismarkAln(
        list.files,
	assembly='hg18',
	sample.id = list('ev1', 'set1', 'by4741'),
	treatment = c(0, 1, 2),
	read.context = 'CpG')

pdf('data/output/methylationStats.pdf')

for (i in 1:length(list.files)){
    getMethylationStats(objs[[i]], plot=T, both.strands=F)}

dev.off()

#meta gene plot
library('GenomicFeatures')

sample.ids = getSampleID(objs)

objs_rg_list = list()
for (i in 1:length(objs)){
   objs_rg_list[[i]] = as(objs[[i]], 'GRanges')
   objs_rg_list[[i]]$score = objs_rg_list[[i]]$numCs/objs_rg_list[[i]]$coverage
}# convert methylRawList object to list of GRange objects

db = makeTxDbFromGFF('data/ref/annotation.gtf', format = 'gtf')#read ref genome annotation and extract gene body; upper and lower flank regions for plotting
gene_rg = genes(db)
gene_rg_upstream = flank(resize(gene_rg, width=1, fix='start'), width=2000)
gene_rg_downstream = flank(resize(gene_rg, width=1, fix='end'), width=2000, start=F)

#plot for gene body
mx_gene = ScoreMatrixList(objs_rg_list, gene_rg, weight.col='score', bin.num=20, strand.aware=T, is.noCovNA=T)
sub_mx_gene = intersectScoreMatrixList(mx_gene, reorder=F)

pdf('data/output/metaGenePlot_gene_body.pdf')
plotMeta(sub_mx_gene, line.col = rainbow(4), profile.names=getSampleID(objs), main='metaGene plot for gene body')

#plot for upstream flank
mx_Up = ScoreMatrixList(objs_rg_list, gene_rg_upstream, weight.col='score', bin.num=20, strand.aware=T, is.noCovNA=T)
sub_mx_Up = intersectScoreMatrixList(mx_Up, reorder=F)

pdf('data/output/metaGenePlot_gene_upstream.pdf')
plotMeta(sub_mx_Up, line.col = rainbow(4), profile.names=getSampleID(objs), main='metaGene plot for upstream of gene')

#plot for downstream flank
mx_Down = ScoreMatrixList(objs_rg_list, gene_rg_downstream, weight.col='score', bin.num=20, strand.aware=T, is.noCovNA=T)
sub_mx_Down = intersectScoreMatrixList(mx_Down, reorder=F)

pdf('data/output/metaGenePlot_gene_downstream.pdf')
plotMeta(sub_mx_Down, line.col = rainbow(4), profile.names=getSampleID(objs), main='metaGene plot for downstream of gene')
