output_path = 'data/trimmed/{sample}_trimmed.fq'


rule trim_galore:
    input:
        'data/{sample}/{sample}.fastq'
    output:
        output_path
    conda:
        '../env/trim_galore.yaml'
    params:
        length = config['trimmed_length'],
        output_dir = os.path.dirname(output_path)
    shell:
        'trim_galore --length {params.length} {input} -o {params.output_dir}'
