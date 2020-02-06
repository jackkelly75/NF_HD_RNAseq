#!/usr/bin/env nextflow

/*
 * Defines pipeline parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */
params.outdir = "results"
params.bins = "/home/jack/Downloads/temp/Qreports/*"

log.info """\
 R N A S E Q - N F   P I P E L I N E
 ===================================
 outdir       : ${params.outdir}
 bins         : ${params.bins}
 """


process Sreports {
    tag "$Sreports"
    publishDir "/home/jack/Downloads/temp/Sreports"

    input:

    output:
    file "*.html" into htmlrecordsPre

    script:
    """
    Sreport -i /home/jack/Downloads/temp/Qreports/ -t Q -o my_test_summary_report
    """
}


process SreportsPost {
    tag "$SreportsPost"
    publishDir "/home/jack/Downloads/temp/Sreports"

    input:

    output:
    file "*.html" into htmlrecordsPost

    script:
    """
    Sreport -i /home/jack/Downloads/temp/reports/ -t P -o results_after_filtering
    """
}
