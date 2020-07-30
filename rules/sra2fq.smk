rule sra2fq:
    input:
        'data/{sample}.sra'
    output:
        'data/{sample}/{sample}.fastq'
    conda:
        '../env/fastq-dump.yaml'
    shell:
        'fastq-dump -O data/{wildcards.sample} {input}'

