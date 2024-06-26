process.executor = 'slurm'
executor.$local.cpus = 1
executor.$local.memory = '32 GB'
process.cpus = 32
process.queue = 'cluster'
process.memory = '64 GB'

params.outdir='/home/edina/testing'

process {
    withName:fastqc {
        container = 'file:///opt/singularity_images/fastqc:0.11.9--hdfd78af_1.sif'
    }
   
}
singularity {
    enabled = true
}


params {
 read = '/home/edina/testing/SRR098034.fastq'
outdir = '/home/edina/testing
}