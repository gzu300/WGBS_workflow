rule bismark_convert:
	input:
		'data/sorted_index/{sample}_sorted.bam'
	output:
		'data/bismark_report/{sample}/{sample}_sorted.bismark.cov.gz'
	conda:
		'../env/bismark.yaml'
	params:
		output_dir = 'data/bismark_report',
		ref_dir = 'data/ref'
	shell:
		'bismark_methylation_extractor --bedGraph --cytosine_report -o {params.output_dir} --genome_folder {params.ref_dir} -s {input}'

rule coverage2cytosine:
	input:
		'data/bismark_report/{sample}/{sample}_sorted.bismark.cov.gz'
	output:
		'data/bismark_report/meth_info_{sample}.CpG_report.txt'
	conda:
		'../env/bismark.yaml'
	params:
		ref_dir = 'data/ref'
	shell:
		'coverage2cytosine -o {output} --genome_folder {params.ref_dir} {input}'
