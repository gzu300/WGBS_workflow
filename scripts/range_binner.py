#!/usr/bin/env python 
'''
put meth coordinates into bins.
'''
import pdb
import pandas as pd
import statistics as stats
import argparse

def parse_gtf_to_dict(gtf_iter):
    #print(each_line.split())
    chrom, _, biotype, start, end, _, strand, _, _, gene_id, *_ = gtf_iter.split()
    return dict(zip(['chrom', 'start', 'end', 'strand', 'gene_id', 'biotype'], 
                    [chrom, int(start), int(end), strand, gene_id[1:-2], biotype]))

def parse_coord_to_dict(coord_iter):
    chrom, coord, strand, numCs, NotCs, *_ = coord_iter.split()
    numCs, NotCs = float(numCs), int(NotCs)
    score = numCs/(numCs+NotCs) if (numCs+NotCs) != 0 else None
    return dict(zip(['chrom', 'coord', 'strand', 'score'], 
                    [chrom, int(coord), strand, score]))

def bin_coords(gene, meth_info, bin_num):
    gene_size = gene['end'] - gene['start']
    bin_size = gene_size//bin_num
    bin_size = bin_size if bin_size > 0 else 1
    if meth_info['strand'] == '+':
        nth_bin = (meth_info['coord'] - gene['start'])//bin_size  #0 based index
                                                                  #so no need to plus 1
    
    else:
        nth_bin = (gene['end'] - meth_info['coord'])//bin_size - 1
    if nth_bin >= bin_num: return bin_num - 1
    if nth_bin < 0: return 0
    return nth_bin

def bin_worker(gtf_file, coord_file, bin_num):
    with open(gtf_file, 'r') as gtf:
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
            with open(coord_file, 'r') as meth:
                score_list = []
                for _ in range(bin_num):
                    score_list.append([])
                for each_coord in meth:
                    coord = parse_coord_to_dict(each_coord)
                    if coord['score'] == None: #this includes score is 0 and None
                                           #0 and None are all treated as 0 in
                                           #the end
                        continue
                    if coord['strand'] == gene['strand']:
                        if coord['coord'] < gene['start']: continue
                        if coord['coord'] > gene['end']: break
                        nth_bin = bin_coords(gene, coord, bin_num)
                        score_list[nth_bin].append(coord['score'])
                yield (gene['gene_id'], score_list)
def calc_ave_score(bin_dict):
    for gene_id, bin_list in bin_dict:
        ave_score_list = []
        for score_list in bin_list:
            if not score_list:
                ave_score_list.append(0)
                continue
            ave_score_list.append(stats.mean(score_list))
        yield (gene_id, ave_score_list)

def main():
    args = argparse.ArgumentParser()
    args.add_argument('-g', '--gtf')
    args.add_argument('-c', '--coord')
    args.add_argument('-b', '--bin_num', type=int)
    args.add_argument('-o', '--output')
    parser = args.parse_args()
    bin_dict = bin_worker(parser.gtf, parser.coord, parser.bin_num)
    score_tuple = calc_ave_score(bin_dict)
    df = pd.DataFrame({
        gene_id: score_list for gene_id, score_list in score_tuple
        })
    df.to_csv(parser.output, header=True, sep='\t')
if __name__ == '__main__': main()


