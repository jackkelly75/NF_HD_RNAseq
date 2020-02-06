#!/usr/bin/env nextflow

/*
 * Defines pipeline parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */
params.reads = "/home/jack/Downloads/temp/reports/*{1,2}_good.fq.gz"
params.transcriptome = "/home/jack/Downloads/temp/hsapien.fa.gz"
params.outdir = "results"

log.info """\
 R N A S E Q - N F   P I P E L I N E
 ===================================
 transcriptome: ${params.transcriptome}
 reads        : ${params.reads}
 outdir       : ${params.outdir}
 """


 /*
 * The reference transcriptome file
 */
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
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}  



process quant {
    tag "$pair_id"
    publishDir '/home/jack/Downloads/temp/data'

    input:
    file transcriptome from transcriptome_file
    file index from transcriptome_index
    set pair_id, file(reads) from read_pairs_ch

    output:
    file(pair_id) into quant_ch

    script:
    """
    salmon quant -l A --threads $task.cpus -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id --numBootstraps 30 --validateMappings --maxMMPExtension 7 --seqBias --gcBias
    """
}
