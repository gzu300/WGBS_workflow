output_path = 'data/deduplicated/{sample}.deduplicated.bam'

rule deduplicate:
    input:
        'data/bismark/{sample}_trimmed_bismark_bt2.bam'
    output:
        output_path
    conda:
        '../env/bismark.yaml'
    params:
        output_dir = os.path.dirname(output_path)
    shell:
        'deduplicate_bismark -s -o {wildcards.sample} --output_dir {params.output_dir} --bam {input}'
