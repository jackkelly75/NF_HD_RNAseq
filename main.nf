#!/usr/bin/env nextflow

params.transcriptome = "$baseDir/data/hsapien.fa.gz"
params.reads = "$baseDir/data/*_{1,2}.fastq.gz"
params.outdir = "results"


log.info """\
 N F - H D - R N A S E Q  P I P E L I N E
 ===================================
 reads        : ${params.reads}
 outdir       : ${params.outdir}
 """

transcriptome_file = file(params.transcriptome)

Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { read_pairs_ch; read_pairs2_ch }


process buildIndex {
    tag "$transcriptome.simpleName"

    input:
    file transcriptome from transcriptome_file

    output:
    file 'index' into transcriptome_index

    script:
    """
    salmon index -t $transcriptome -i index -k 31
    """
}  

process trimFilter {
    
    tag "$trimFilter"
    publishDir "1_FastQPuri"

    input:
    tuple val(pair_id), path(reads) from read_pairs_ch

    output:
    file('*good.fq.gz') into goodfiles
    
    script:
    """
    trimFilterPE -f ${reads[0]}:${reads[1]}  -l 101 --trimQ ENDSFRAC --trimN ENDS -m 31 -o $pair_id
    """
}
