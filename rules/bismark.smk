rule build:
    input:
        'data/ref/'
    output:
        'data/ref/Bisulfite_Genome'
    conda:
        '../env/bismark.yaml'
    shell:
        'bismark_genome_preparation {input}'

rule bismark:
    input:
        fq = 'data/trimmed/{sample}_trimmed.fq',
        genome_folder = 'data/ref/'
    output:
        'data/bismark/{sample}_trimmed_bismark_bt2.bam'
    conda:
        '../env/bismark.yaml'
    log: '../log/{sample}_bismark.log'
    shell:
        'bismark {input.genome_folder} {input.fq} -o data/bismark/ 2>{log}'
