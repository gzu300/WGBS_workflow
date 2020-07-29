rule download_reference:
    output:
        genome = 'data/ref/ref_genome.fa.gz'
    shell:
        'wget -O {output.genome} ftp://ftp.ensembl.org/pub/release-{config[genome_release]}/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_rm.toplevel.fa.gz'

rule download_annotation:
    output:
        annotation = 'data/ref/annotation.gtf.gz'
    shell:
        'wget -O {output.annotation} ftp://ftp.ensembl.org/pub/release-{config[genome_release]}/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.{config[genome_release]}.gtf.gz'

rule unzip:
    input:
        'data/ref/ref_genome.fa.gz',
        'data/ref/annotation.gtf.gz'
    output:
        'data/ref/ref_genome.fa',
        'data/ref/annotation.gtf'
    shell:
        'gunzip {input[0]};'
        'gunzip {input[1]}'
