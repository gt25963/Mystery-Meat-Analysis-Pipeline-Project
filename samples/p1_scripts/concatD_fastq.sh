#!/bin/bash

#Setup
##Define which sample I am working with. Allows for me to easily reuse code for other samples. 
sample_name="sampleD"

##Define the output file, following FASTA format
output="${sample_name}_complete.fasta"

##Remove any existing old output file if it exists, to prevent appending old data
rm -f "$output"

#Show progess
##Print where I am, to check if in the correct directory
echo "Current directory: $(pwd)"
echo ""
echo "Concatenating FASTQ files for ${sample_name}"

#Concatenate Sequences Together
##Set up variable to store the combined sequence for sampleA. Parts need to be appended because they have been split into 3 files.
full_sequence=""

##Loop through each of the 3 FASTQ files
##FASTQ files contain the follwoing 4 lines: Name, DNA, +, Quality
for file in ${sample_name}_part1.FASTQ ${sample_name}_part2.FASTQ ${sample_name}_part3.FASTQ
do
    if [ -f "$file" ]; then #Check if the file exists before reading it
        echo "Processing: $file"

        ##Pull out just the DNA sequence from the FASTQ file
        ##NR%4==2 means 'retrieve every 4th line starting from line 2'
        ##tr -d removes spaces and line breaks 
        sequence=$(awk 'NR%4==2' "$file" | tr -d 'n' | tr -d ' ')

        ##Append to full sequence
        #Rrebuilds the complete DNA from the 3 split sequences
        full_sequence="${full_sequence}${sequence}"

        ##Count how many bases have been added
        ## ${variable} indicated how long the sequence is
        base_count=${#sequence}
        echo "  Added${base_count} bases to sequence"
    fi
done

#Save the Result
##Count total bases in concatenated sequence (check all data is present)
total_bases=${#full_sequence}
echo ""
echo "${sample_name}: ${total_bases} bases"

##Save to FASTA format 
## > indicates the name line (header)
echo ">${sample_name}" > "$output"
## >> Puts the DNA sequence on the next line
echo "$full_sequence" >> "$output"

echo "Sequence written to $output"
