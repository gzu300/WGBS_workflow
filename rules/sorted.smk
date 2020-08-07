rule sort:
    input:
        'data/deduplicated/{sample}.deduplicated.bam'
    output:
        'data/sorted_index/{sample}_sorted.bam',
        'data/sorted_index/{sample}_sorted.bam.bai'
    conda: '../env/samtools.yaml'
    shell:
        'samtools sort -o {output[0]} {input} && '
        'samtools index {output[0]}'
