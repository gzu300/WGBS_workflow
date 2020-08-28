#!/usr/bin/env python 
import pdb
import pandas as pd
import statistics as stats

'''
put meth coordinates into bins.
'''
BIN_NUM = 20

def parse_gtf_to_dict(gtf_iter):
    #print(each_line.split())
    chrom, _, biotype, start, end, _, strand, _, _, gene_id, *_ = gtf_iter.split()
    return dict(zip(['chrom', 'start', 'end', 'strand', 'gene_id', 'biotype'], 
                    [chrom, int(start), int(end), strand, gene_id[1:-2], biotype]))

def parse_coord_to_dict(coord_iter):
    chrom, coord, strand, numCs, coverage, *_ = coord_iter.split()
    numCs, coverage = float(numCs), int(coverage)
    score = numCs/coverage if coverage != 0 else None
    return dict(zip(['chrom', 'coord', 'strand', 'score'], 
                    [chrom, int(coord), strand, score]))

def bin_coords(gene, meth_info, bin_num = BIN_NUM):
    gene_size = gene['end'] - gene['start']
    bin_size = gene_size//BIN_NUM
    bin_size = bin_size if bin_size > 0 else 1
    #if gene['strand'] == meth_info['strand']:
#        if (meth_info['coord'] > gene['start'])& \
#           (meth_info['coord'] < gene['end']):
    if meth_info['strand'] == '+':
        nth_bin = (meth_info['coord'] - gene['start'])//bin_size  #0 based index
        # so no need to plus 1
        return nth_bin if nth_bin < BIN_NUM else BIN_NUM-1
    nth_bin = (gene['end'] - meth_info['coord'])//bin_size - 1
    if nth_bin >= BIN_NUM: return BIN_NUM - 1
    if nth_bin < 0: return 0
    return nth_bin

with open('example/anno_chr1.gtf', 'r') as gtf:
    '''
    parse gtf to {'gene_id':{'chrom':xxx, ...}}
    parse meth to {'chrom':xxx, 'coord':xxx, 'score':xxx}
    bin meth score to corresponding coord[[],[],[]...]
    after looping through that gene, calculate average meth score.
    [float, float, ...]
    add to meta collector
    the final output is a dictionary. {'gene_id': [bin1, bin2...bin20]}
    '''
    meta_collector = {}
    for each_gene in gtf:
        gene = parse_gtf_to_dict(each_gene)
        if gene['biotype'] != 'gene': continue
        score_list = []
        for _ in range(BIN_NUM):
            score_list.append([])
        with open('example/meth_info_chrI.txt', 'r') as meth:
            for each_coord in meth:
                coord = parse_coord_to_dict(each_coord)
                #print(coord['strand'], coord['coord'], gene['start'], gene['end'])
                if (coord['coord'] > gene['end']) or \
                   (coord['coord'] < gene['start']) or \
                   (coord['strand'] != gene['strand']):
                    continue

                nth_bin= bin_coords(gene, coord)
#                if nth_bin:
#                    print(nth_bin, coord['strand'], coord['coord'], gene['start'], gene['end'])
                if nth_bin:
                    score_list[nth_bin].append(coord['score'])
            new_list = []
            for each_bin in score_list:
                try:
                    new_list.append(stats.mean(each_bin))
                except (stats.StatisticsError, TypeError):
                    new_list.append(None)
            meta_collector[gene['gene_id']] = new_list
    df = pd.DataFrame(meta_collector, columns=meta_collector.keys())
    print(df.dropna(axis=1, how='all'))


