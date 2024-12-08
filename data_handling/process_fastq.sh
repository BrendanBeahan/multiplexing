#!/bin/bash

# Load the required modules
module load minimap2/2.26-GCCcore-12.3.0
module load SAMtools/1.18-GCC-12.3.0

# Define the reference FASTA file
REFERENCE="synthetic_sequences.fasta"

# Set the number of threads
THREADS=8

# Initialize the read counts output file
READ_COUNTS_FILE="read_counts.txt"
echo -e "Barcode\tReference\tRead_Counts" > $READ_COUNTS_FILE

# Loop through each FASTQ file in the directory
for INPUT_FASTQ in barcode_*.fastq; do
    # Extract the number from the filename
    NUMBER=$(echo $INPUT_FASTQ | grep -oP '\d+(?=\.fastq)')
    
    echo "Processing $INPUT_FASTQ (number $NUMBER)..."
    
    # Step 1: Align with minimap2
    SAM_FILE="alignments_$NUMBER.sam"
    LOG_FILE="minimap2_$NUMBER.log"
    minimap2 -x map-ont -t $THREADS -a $REFERENCE $INPUT_FASTQ > $SAM_FILE 2> $LOG_FILE
    if [[ ! -s $SAM_FILE ]]; then
        echo "Error: Failed to generate $SAM_FILE for $INPUT_FASTQ. Check $LOG_FILE for details. Skipping this file."
        continue
    fi

    # Step 2: Convert SAM to BAM, sort, and save
    BAM_FILE="alignments_$NUMBER.sorted.bam"
    samtools view -bS $SAM_FILE | samtools sort -@ $THREADS -o $BAM_FILE
    if [[ ! -s $BAM_FILE ]]; then
        echo "Error: Failed to generate $BAM_FILE for $INPUT_FASTQ. Skipping this file."
        continue
    fi

    # Step 3: Filter for primary alignments
    PRIMARY_SAM_FILE="primary_alignments_$NUMBER.sam"
    samtools view -F 2304 -h $BAM_FILE > $PRIMARY_SAM_FILE
    if [[ ! -s $PRIMARY_SAM_FILE ]]; then
        echo "Error: Failed to generate $PRIMARY_SAM_FILE for $INPUT_FASTQ. Skipping this file."
        continue
    fi

    # Step 4: Convert filtered SAM to BAM, sort, and index
    PRIMARY_BAM_FILE="primary_alignments_$NUMBER.bam"
    samtools view -bS $PRIMARY_SAM_FILE | samtools sort -o $PRIMARY_BAM_FILE
    samtools index $PRIMARY_BAM_FILE
    if [[ ! -s $PRIMARY_BAM_FILE ]]; then
        echo "Error: Failed to generate $PRIMARY_BAM_FILE for $INPUT_FASTQ. Skipping this file."
        continue
    fi

    # Step 5: Generate idxstats
    IDXSTATS_FILE="primary_alignments_$NUMBER.idxstats"
    samtools idxstats $PRIMARY_BAM_FILE > $IDXSTATS_FILE
    if [[ ! -s $IDXSTATS_FILE ]]; then
        echo "Error: Failed to generate $IDXSTATS_FILE for $INPUT_FASTQ. Skipping this file."
        continue
    fi

    # Parse idxstats to append to read_counts.txt
    echo "Parsing idxstats for $INPUT_FASTQ..."
    while read -r LINE; do
        REF_NAME=$(echo "$LINE" | awk '{print $1}')
        READ_COUNT=$(echo "$LINE" | awk '{print $3}')
        if [[ $REF_NAME == "*" ]]; then
            REF_NAME="unmapped"
        fi
        echo -e "$NUMBER\t$REF_NAME\t$READ_COUNT" >> $READ_COUNTS_FILE
    done < $IDXSTATS_FILE

    echo "Finished processing $INPUT_FASTQ."
done

echo "All files processed."
