#!/bin/bash
#Define the sample name
sample_name="sampleD"

#Define the output file
output="${sample_name}_complete.fasta"

#Remove old output file if it exists
rm -f "$output"

echo "Current directory: $(pwd)"
echo ""
echo "Concatenating FASTQ files for ${sample_name}"

#Set up variable to store sampleA
full_sequence=""

#Loop through the 3 FASTQ files
for file in ${sample_name}_part1.FASTQ ${sample_name}_part2.FASTQ ${sample_name}_part3.FASTQ
do
    if [ -f "$file" ]; then
        echo "Processing: $file"

        #Extract sequences from FASTQ file, every 4th line starting from line 2
        sequence=$(awk 'NR%4==2' "$file" | tr -d 'n' | tr -d ' ')

        #Append to full sequence
        full_sequence="${full_sequence}${sequence}"

        #Count bases for this part of sample
        base_count=${#sequence}
        echo "  Added${base_count} bases to sequence"
    fi
done

#Count total bases in concatenated sequence
total_bases=${#full_sequence}
echo ""
echo "${sample_name}: ${total_bases} bases"

#Write the header and sequence to output
echo ">${sample_name}" > "$output"
echo "$full_sequence" >> "$output"

echo "Sequence written to $output"
