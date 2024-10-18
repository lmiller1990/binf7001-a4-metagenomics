#!/bin/bash

prefixes=("Bran" "CO" "PS")
root="/opt/BINF7001/2023/Prac11_2023/"
genome_index="$root/index_files/GCA_014282415.2.fa"


# Loop over each prefix
for prefix in "${prefixes[@]}"; do
  # Loop over each number from 1 to 3
  for j in {1..3}; do
    # Define file names for forward and reverse reads
    forward="$root/fastq/${prefix}_${j}_R1.truncated.fastq.gz"
    reverse="$root/fastq/${prefix}_${j}_R2.truncated.fastq.gz"

    # Define output file names for non-host reads
    output_forward="${prefix}_${j}_R1_nonhost.fastq"
    output_reverse="${prefix}_${j}_R2_nonhost.fastq"

    # Align reads with bwa mem and process with samtools and bamToFastq
    bwa mem -t 4 "$genome_index" "$forward" "$reverse" | \
    samtools view -b -f 12 -F 256 -F 2048 | \
    bamToFastq -i /dev/stdin -fq "$output_forward" -fq2 "$output_reverse"
  done
done
