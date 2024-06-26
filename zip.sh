#!/bin/sh

echo $1
file=$1

if [[ $file = *.zip ]]
then
    echo "Unzipping: $file"
	sample="$(basename $1 .zip)"
	mkdir ./$sample
    unzip $file -d ./$sample
    relative_sample_path="./$sample"
    absolute_sample_path="$(cd "$(dirname $relative_sample_path)"; pwd)/$(basename $relative_sample_path)"
    files="$absolute_sample_path/${sample}_{1,2}.fastq"
	echo $files
    module load bioconda
    nextflow -C /home/edina/testdirectory/runs/check.config run /home/edina/testdirectory/runs/raw_data_analysis.nf --reads $files
else
    echo "No zip on input!"
fi