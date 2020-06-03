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


process trimFilter {
    
    tag "$trimFilter"
    publishDir "1_FastQPuri"

    input:
    file(reads) from read_pairs_ch

    output:
    file('*{1,2}_good.fq.gz') into goodfiles
    
    script:
    """
    trimFilterPE -f ${reads[0]}:${reads[1]}  -l 101 --trimQ ENDSFRAC --trimN ENDS -m 31 -o ${reads[0].baseName}
    """
}
