#!/usr/bin/bash

# Loop creating job script for every TCGA file

#Args:
#manifest   Directory to the gdc manifest or missing files from step 02
#token      Directory to the gdc user token for controlled data
#bampath    Directory to bam files downloaded from gdc
#pathout    Path to desired output

#$ -e ./$JOB_NAME-$JOB_ID.err
#$ -o ./$JOB_NAME-$JOB_ID.log

manifest="/data/projects/2020/OvarianCancerHH/Thesis_Katja/PNAT_Results/missing_files.txt"
token="/data/projects/2020/OvarianCancerHH/Thesis_Katja/gdc-user-token.2024-04-03T080557.058Z.txt"
bampath="/data/projects/2020/OvarianCancerHH/Thesis_Katja/PNAT_RNAseq/"
pathout="/data/projects/2020/OvarianCancerHH/Thesis_Katja/PNAT_Results/"

# Job count, to submit in batches
max_concurrent_jobs=3

get_running_job_count() {
    qstat -u 'rungger' | grep 'rungger' | wc -l
}

# Loop over manifest file, removing the header with tail and only reading in the first column

tail -n +2 "$manifest" | while IFS= read -r line; do

    while [ $(get_running_job_count) -ge $max_concurrent_jobs ]; do
        sleep 5 
    done

    id=$(echo "$line" | awk '{print $1}')
    echo "[+] qsub for ID: $id"
    qsub -q all.q@apollo-01.local /home/rungger/myScratch/identification_bacteria/NEW_Pipe/01b_combi_down_ident.sh $bampath $pathout $id $token

done
