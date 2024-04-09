#!/bin/sh
#$ -S /bin/sh
#####source /etc/profile.d/modules.sh
#$ -pe smp 5
#$ -cwd
#$ -V
#### Jobdescription at qstat
#$ -N Create_BIOM
#### Error Outputfile
#$ -e ./$JOB_NAME-$JOB_ID.err
#$ -o ./$JOB_NAME-$JOB_ID.log
#### Resubmit
#$ -r y

input="/data/projects/2020/OvarianCancerHH/Thesis_Katja/PNAT_Results"
bampath="/data/projects/2020/OvarianCancerHH/Thesis_Katja/PNAT_RNAseq"

cd $bampath || exit
count=$(find . -mindepth 1 -maxdepth 1 -type d | wc -l)
echo "Number of directories: $count"


source activate KrakenTools
python3     /home/rungger/myScratch/identification_bacteria/NEW_Pipe/03_create_biom.py --input $input --n_samples $count

