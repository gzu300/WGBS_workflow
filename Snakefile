configfile: 'config.yaml'

rule all:
    input:
        expand('data/bismark/{sample}_trimmed_bismark_bt2.bam', sample=config['samples']),
        'data/ref/ref_genome.fa',
        'data/ref/annotation.gtf'

rule download:
    output:
        'data/{sample}.sra'
    params: lambda wildcards: config['samples'][wildcards.sample]
    shell:
        'wget -O {output} https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/{params}/{params}.1' 

#rule sra2fq:
#    input:
#        'data/{sample}.sra'
#    output:
#        'data/{sample}/{sample}.fastq'
#    conda:
#        'env/fastq-dump.yaml'
#    shell:
#        'fastq-dump -O data/{wildcards.sample} {input}'

include: 'rules/sra2fq.smk'
include: 'rules/download_ref.smk'
include: 'rules/trim_galore.smk'
include: 'rules/bismark.smk'
