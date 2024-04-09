#!/bin/sh
#$ -S /bin/sh
#####source /etc/profile.d/modules.sh
#$ -pe smp 5
#$ -cwd
#$ -V
#### Jobdescription at qstat
#$ -N Extract_reads
#### Error Outputfile
#$ -e ./$JOB_NAME-$JOB_ID.err
#$ -o ./$JOB_NAME-$JOB_ID.log
#### Resubmit
#$ -r y

species="2719117"
barcode_path="/data/projects/2020/OvarianCancerHH/Thesis_Katja/OV_RNAseq/TCGA_barcodes.txt"
result_path="/data/projects/2020/OvarianCancerHH/Thesis_Katja/OV_Results_NEW"
pathout="/data/projects/2020/OvarianCancerHH/Thesis_Katja/OV_Data_nt/Genome_alignment"

# Create new file to append read number
cd $pathout || exit
name_output=$species"_reads.txt"
name_output_path=$pathout/$name_output

echo "Barcode Shortened UUID Read_number" > $name_output 


# Loop over barcode file, removing the header with tail and only reading in the first column
echo "[+] Starting loop to get sequences"


tail -n +2 "$barcode_path" | while IFS= read -r line; do

    barcode=$(echo "$line" | awk '{print $1}')
    barcode_short=${barcode:0:16} 

    kraken_file=$barcode_short".kraken.txt"
    kraken_path=$result_path/$kraken_file

    # Extracting read numbers from kraken.txt output
    read_number=$(awk -v species="$species" '$3 == species { print $2}' "$kraken_path")

    # Temporary file for read ID of one sample
    temp_file=$pathout/$barcode_short"_reads.txt"

    # Check if there are any reads
    if [ -z "$read_number" ] 
    then
        echo "- No reads for $species in $barcode_short"

    # Outputting sample and read information
    # Writing read id to secondary temporary file to be faster in samtools view
    else
        echo "+ Getting sequences for tax id $species for file $kraken_file"
        for read in $read_number
        do
        echo "${line} ${read}" >> $name_output_path
        echo "${read}" >> "$temp_file"
        done 

        # Creating txt file to output sequences
        cd $pathout || exit
        name_seq=$species"_sequences.bam"
        name_seq_path=$pathout/$name_seq

        # Path to bam file with unaligned reads
        bam_file="${barcode_short}_merged.bam"
        bam_path="$result_path/$bam_file" 

        # Extract the reads
        samtools view "$bam_path" | awk 'BEGIN {while(getline < "'"$temp_file"'") reads[$1]=1} $1 in reads {print $0}' >> $name_seq_path

        # Remove temp_file
        rm $temp_file
    
    fi

done

# Create a minimal header file
header_sam=$pathout"/header.sam"
echo -e "@HD\tVN:1.6\tSO:unsorted" > $header_sam 
echo -e "@SQ\tSN:reference\tLN:1000" >> $header_sam  

# Combine the header and your BAM file
name_seq=$species"_sequences.bam"
name_seq_path=$pathout/$name_seq
new_bam=$pathout/$species"_header.bam"
cat $header_sam $name_seq_path > $new_bam

## Converting the bam to fasta
fasta_output=$species"_allreads.fasta"
fasta_output_path=$pathout/$fasta_output

samtools fasta $new_bam > $fasta_output_path 

rm $header_sam

