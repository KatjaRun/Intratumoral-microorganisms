#!/bin/sh
#$ -S /bin/sh
#####source /etc/profile.d/modules.sh
#$ -pe smp 5
#$ -cwd
#$ -V
#### Jobdescription at qstat
#$ -N STAR_bacteria
#### Error Outputfile
#$ -e ./$JOB_NAME-$JOB_ID.err
#$ -o ./$JOB_NAME-$JOB_ID.log
#### Resubmit
#$ -r y

species="2719117"
pathout="/data/projects/2020/OvarianCancerHH/Thesis_Katja/OV_Data_nt/Genome_alignment"
path_download=$pathout/$species"_dataset"

# Specifying paths for running STAR
genome_dir=$path_download"/STAR/"
genome_fasta=$path_download"/cds_from_genomic.fna"
genome_gtf=$path_download"/genomic.gtf"
fasta=$pathout/$species"_allreads.fasta"
prefix=$pathout/$species"_STAR"
counts=$pathout/$species"_counts.txt"

## Creating a STAR genome index
STAR --runMode genomeGenerate --runThreadN 5 --genomeDir "$genome_dir" --genomeFastaFiles "$genome_fasta" #--sjdbGTFfile "$genome_gtf"

## Run the mapping process
STAR --runThreadN 5 --genomeDir "$genome_dir" --readFilesIn "$fasta" --outFileNamePrefix "$prefix" --genomeSAindexNbases 9 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0

## Running feature counts
star_out=$prefix"Aligned.out.sam"

echo $star_out
featureCounts -a "$genome_gtf" -o $counts $star_out

## Sorting samtools output
samtools sort $star_out -O sam -o $star_out

## Calculating depth
depth=$species"_CDS_samtools.coverage"
depth_path=$pathout/$depth

samtools depth -a $star_out > $depth_path

