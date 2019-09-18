#!/usr/bin/env nextflow

/*
 * Defines pipeline parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */
params.reads = "/home/jack/Downloads/temp/*_{1,2}.fastq.gz"
params.outdir = "results"


log.info """\
 R N A S E Q - N F   P I P E L I N E
 ===================================
 reads        : ${params.reads}
 outdir       : ${params.outdir}
 """


Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { read_pairs_ch; read_pairs2_ch }




Channel
    .fromPath('file_info_48.csv')
    .splitCsv()
    .set { singleReads48 }

Channel
    .fromPath('file_info_48_name.csv')
    .splitCsv()
    .set { singleReads48name }

Channel
    .fromPath('file_info_96.csv')
    .splitCsv()
    .set { singleReads96 }

Channel
    .fromPath('file_info_96_name.csv')
    .splitCsv()
    .set { singleReads96name }


process Qreports48 {
    tag "$Qreport48"
    publishDir "/home/jack/Downloads/temp/Qreports"

    input:
    set name, file(reads) from singleReads48
    set names, file(eachname) from singleReads48name

    output:
    file "*.bin" into binfiles48
    file "*.html" into htmlfiles48

    script:
    """
    Qreport -i $name -l 101 -t 48 -o $names
    """
}




process Qreports96 {
    tag "$Qreport96"
    publishDir "/home/jack/Downloads/temp/Qreports"

    input:
    set name, file(reads) from singleReads96
    set names, file(eachname) from singleReads96name

    output:
    file "*.bin" into binfiles96
    file "*.html" into htmlfiles96

    script:
    """
    Qreport -i $name -l 101 -t 96 -o $names
    """
}





process trimFilter {
    tag "$trimFilter"
    publishDir "/home/jack/Downloads/temp/reports"

    input:
    set pair_id, file(reads) from read_pairs_ch


    output:
    file "*.bin" into fastqbinfiles
    file "*good.fq.gz" into goodfiles

    script:
    """
    trimFilterPE -f ${reads[0]}:${reads[1]}  -l 101 --trimQ ENDSFRAC --trimN ENDS -m 25 -o $pair_id
    """
}