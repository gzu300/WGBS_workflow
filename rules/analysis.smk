#print(expand('data/sorted_index/{sample}_sorted.bam', sample=config['samples'])


rule stats:
    input: 
        expand('data/sorted_index/{sample}_sorted.bam', sample=config['samples'])

    output:
        'data/output/methylationStats.pdf'
    params:
        config['samples']
    conda:
        '../env/r-base.yaml'
    log:
        'log/stats.log'
    shell:
        'Rscript scripts/stats.R'

