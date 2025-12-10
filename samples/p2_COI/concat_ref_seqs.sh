#!/bin/bash
#Combine all 16 reference sequences into one long fasta file
cat *.fasta > FINAL_all_ref_sequences.fasta

#Error Handling: check if file was created and is not empty
#[ ! -s filename ] = checks if file doesn't exist OR has zero size
#echo "message" >&2 = prints to standard error
#exit 1 = exits with error code
if [ ! -s FINAL_all_ref_sequences.fasta ]; then
    echo "ERROR: No sequences concatenated or file is empty!" >&2
    exit 1
fi
