#!/usr/bin/env nextflow

if (!params.fastq){params.fastq = ""} 
read_ch = file(params.fastq, type: 'any')

if (!params.outdir){params.outdir = ""} 
out_ch = file(params.outdir, type: 'any')

process fastqc {

input:
 file read from read_ch
 
output:
 file "./${read.simpleName}_results/${read.simpleName}_fastqc/summary.txt" into out_ch


"""
mkdir ./${read.simpleName}_results
fastqc -t 16 -o ./${read.simpleName}_results --extract \$PWD/${read}
"""
}



workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}

