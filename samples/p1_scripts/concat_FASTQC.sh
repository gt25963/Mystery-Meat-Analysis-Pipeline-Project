#!/bin/bash

echo "Concatenating FASTQ files..."

cat sampleA_part1.FASTQ sampleA_part2.FASTQ sampleA_part3.FASTQ > sampleA.fastq
echo "✓ sampleA.fastq created"

cat sampleB_part1.FASTQ sampleB_part2.FASTQ sampleB_part3.FASTQ > sampleB.fastq
echo "✓ sampleB.fastq created"

cat sampleC_part1.FASTQ sampleC_part2.FASTQ sampleC_part3.FASTQ > sampleC.fastq
echo "✓ sampleC.fastq created"

cat sampleD_part1.FASTQ sampleD_part2.FASTQ sampleD_part3.FASTQ > sampleD.fastq
echo "✓ sampleD.fastq created"

echo "All files concatenated successfully!"