#!/bin/sh
#$ -S /bin/sh
#####source /etc/profile.d/modules.sh
#$ -cwd
#$ -V
#### Jobdescription at qstat
#$ -N Download_Identification
#### Error Outputfile
#$ -e ./$JOB_NAME-$JOB_ID.err
#$ -o ./$JOB_NAME-$JOB_ID.log
#### Resubmit
#$ -r y

bampath=$1
pathout=$2
id=$3
token=$4

echo "[+] The id being downloaded: $id"

# Download the file 
gdc-client download $id -t $token -d $bampath

# Getting input
input="$bampath/$id/*.bam" 
echo "[+] The input being worked on: $input"

## Script to get barcode out of R
R_SCRIPT='library(TCGAutils);
new <- UUIDhistory("'$id'");
bar <- UUIDtoBarcode(new[which.max(new$data_release), "uuid"], from_type = "file_id");
cat(bar[1, 2])'

# Running R script and extracting the barcode
R_OUTPUT=$(Rscript -e "$R_SCRIPT")
barcode=${R_OUTPUT:0:16}  

echo "[+] The barcode is: $R_OUTPUT, or shortened: $barcode"

# Filtering out unmapped reads
tmpoutput=$pathout/$barcode

echo "[+] The tmp output is: $tmpoutput"

tmps1=$tmpoutput"_tmps1.bam"
tmps2=$tmpoutput"_tmps2.bam"
tmps3=$tmpoutput"_tmps3.bam"

samtools view -u -f 4  -F 264 $input  > $tmps1
samtools view -u -f 8  -F 260 $input  > $tmps2
samtools view -u -f 12 -F 256 $input  > $tmps3

# Merge all tmp files
merged_unaligned=$tmpoutput"_merged.bam"

samtools merge -u - $tmps1 $tmps2 $tmps3 | samtools sort -n - -T tmp_sort -o $merged_unaligned

echo "[+] Created merged unaligned output: $merged_unaligned"

rm $tmpoutput"_tmps1.bam"
rm $tmpoutput"_tmps2.bam"
rm $tmpoutput"_tmps3.bam"

# File names for output
output1=$pathout/$barcode"_R1.fq.gz"
output2=$pathout/$barcode"_R2.fq.gz"

echo "[+] Creating output files: $output1 and $output2 ."

## Creating fastq files
gatk SamToFastq -I $merged_unaligned -F $output1 -F2 $output2 -NON_PF true


## Starting identification pipeline
echo "[+] Starting identification pipeline"

python3  /home/rungger/myScratch/identification_bacteria/NEW_Pipe/01c_pipe_identification.py --path_R1 $output1 --path_R2 $output2 --output $pathout --base $barcode




