#!/usr/bin/env nextflow

if (!params.fastq){params.fastq = ""} 
read_ch = file(params.fastq, type: 'any')

if (!params.fastqs){params.fastqs = ""} 
Channel
  .fromFilePairs( params.reads, checkIfExists: true )
  .set { read_pairs_ch }

process fastqc {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*\/summary.txt$/) "${read.simpleName}/$filename"}

if (mode == 'SE'){

input:
	file read from read_ch

}else if ( mode == 'PE'){

input:
	set val(name), file(reads) from read_pairs_ch
 
output:
	file "./${read}_results/${read.simpleName}_fastqc/summary.txt" into out_ch
file "./${reads}_results/${reads.simpleName}_fastqc/summary.txt"
script:
if ( mode == 'SE')


    """
    mkdir ./${read.simpleName}_results
	fastqc -t 16 -o ./${read.simpleName}_results --extract \$PWD/${read}
    """

    else if ( mode == 'PE')

    """
    mkdir ./${reads.simpleName}_results
	fastqc -t 16 -o ./${reads.simpleName}_results --extract \$PWD/${reads}
    """
workflow{
read_pairs_ch.view()
}

}



workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}

