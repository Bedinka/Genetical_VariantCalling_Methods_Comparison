nextflow.enable.dsl = 2

params.reads = "$baseDir/paired/*_{1,2}.fq"
params.outdir = "results"

log.info """\
 Q U A L I T Y - C H E C K  P I P E L I N E
 ===================================
 reads        : ${params.reads}
 outdir       : ${params.outdir}
 """

include { FASTQ } from './modules/fastqc'

workflow {
  read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true ) 
  FASTQ( params.fastqs, read_pairs_ch
	FASTQ( params.fastq, read_ch )

}

workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}

