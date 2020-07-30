rule trim_galore:
    input:
        'data/{sample}/{sample}.fastq'
    output:
        'data/{sample}_trimmed.fq'
    conda:
        '../env/trim_galore.yaml'
    params:
        length = config['trimmed_length']
    shell:
        'trim_galore --length {params.length} {input} -o data'
