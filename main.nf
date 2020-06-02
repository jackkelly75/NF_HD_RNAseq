#!/usr/bin/env nextflow

params.outdir = "results"
params.files = "$baseDir/data/*_1.fastq.gz"

log.info """\
 N F - H D - R N A S E Q  P I P E L I N E
 ===================================
 outdir       : ${params.outdir}
 files        : ${params.files}	
 """
   
Channel
    .from(params.files)
    .ifEmpty { error "Cannot find any reads matching: ${params.files}" }
    .set { raw_reads_fastqc }
            

process trimFilter {
    
    tag "$trimFilter"
    publishDir "1_FastQPuri"

    input:
    file(reads) from raw_reads_fastqc

    output:
    set val(pair_id), file('*{1,2}_good.fq.gz') into goodfiles
    
    shell:
    '''
	   a=$(echo !{reads} | sed -e 's/_1/_2/')
	   ln=${!{reads}##*/}
   	v2=${ln::-10}
   	trimFilterPE -f !{reads}:$a -l 101 --trimQ ENDSFRAC --trimN ENDS -m 31 -o ${v2}
    '''
}
