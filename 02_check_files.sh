#!/bin/sh
#$ -S /bin/sh
#####source /etc/profile.d/modules.sh
#$ -pe smp 5
#$ -cwd
#$ -V
#### Jobdescription at qstat
#$ -N Check_files
#### Error Outputfile
#$ -e ./$JOB_NAME-$JOB_ID.err
#$ -o ./$JOB_NAME-$JOB_ID.log
#### Resubmit
#$ -r y

r_script="/home/rungger/myScratch/identification_bacteria/NEW_Pipe/02_get_barcodes.R"
bam_path="/data/projects/2020/OvarianCancerHH/Thesis_Katja/PNAT_RNAseq/"
pathout="/data/projects/2020/OvarianCancerHH/Thesis_Katja/PNAT_Results/"

# Run R script
barcode_file=$(Rscript "$r_script" "$bam_path")

echo "[+] Created the file: $barcode_file"

# Create new file to keep original manifest file
cd $pathout || exit
echo "Missing uuids" > missing_files.txt

# Check if translated barcode is present as _bracken_species file, if not append to txt file
tail -n +2 "$barcode_file" | while IFS= read -r line; do

    barcode=$(echo "$line" | awk '{print $2}')
    
    if ! ls "$pathout"/*_bracken_species.kreport2 | grep -q "$barcode"; then
        echo "$line" | awk '{print $3}' >> missing_files.txt
        echo "[+] File of $barcode is missing, adding to txt file"
    fi

done



