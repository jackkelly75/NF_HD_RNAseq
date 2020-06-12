#!/usr/bin/env nextflow

params.transcriptome = "$baseDir/hsapien.fa.gz"
params.reads = "$baseDir/1_FastqPuri/*_{1,2}_good.fq.gz"
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

process quant {
    cpus = 2
    tag "$pair_id"
    publishDir '2_quant'


    input:    
    file index from transcriptome_index
    tuple val(pair_id), path(reads) from read_pairs_ch

    output:
    file(pair_id) into quant_ch

    script:
    """
    salmon quant -l A -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id --validateMappings --seqBias --gcBias
    """
}
