#!/usr/bin/env nextflow

if (!params.reads){params.reads = ""} 

Channel
  .fromFilePairs( params.reads, checkIfExists: true )
  .set { read_pairs_ch }

workflow{
read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true ) 
myFile = file('summary.txt')
myFile.renameTo('${reads.simpleName)_summary.txt')
}

process fastqc {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*\/summary.txt$/) "${read.simpleName}/$filename"}

input:
 set val(name), file(reads) from read_pairs_ch
 
output:
 file "./${reads.simpleName}_results/${reads.simpleName}_fastqc/${reads.simpleName}_summary.txt" into out_ch


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

