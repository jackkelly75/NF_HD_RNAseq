#!/usr/bin/env nextflow

params.transcriptome = "$baseDir/data/hsapien.fa.gz"
params.reads = "$baseDir/data/*_{1,2}.fastq.gz"
params.outdir = "results"
params.files = "$baseDir/data/*_1.fastq.gz"

log.info """\
 N F - H D - R N A S E Q  P I P E L I N E
 ===================================
 reads        : ${params.reads}
 outdir       : ${params.outdir}
 files        : ${params.files}	
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
    tuple val(pair_id), path(reads) from read_pairs_ch

    output:
    set val(pair_id), file('*{1,2}_good.fq.gz') into goodfiles
    
    shell:
    '''
    for fn in !{FILES};
    do
	a=$(echo ${fn} | sed -e 's/_1/_2/')
	ln=${fn##*/}
    	v2=${ln::-10}
    	trimFilterPE -f $fn:$a -l 101 --trimQ ENDSFRAC --trimN ENDS -m 31 -o ${v2}
    done
    '''
}
