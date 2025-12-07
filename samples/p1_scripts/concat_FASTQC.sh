#!/bin/bash
#Join Sample A files 
##Combine the 3 parts of sample A into one complete fastq file
##cat concatenates the files together. rebuilding the sample fastq for FASTQC analysis 
cat sampleA_part1.FASTQ sampleA_part2.FASTQ sampleA_part3.FASTQ > sampleA.fastq
##Confirm the parts have been combined
echo "sampleA.fastq created"

#Join Sample B files
##Concatemate the 3 parts for sample B into one complete file
cat sampleB_part1.FASTQ sampleB_part2.FASTQ sampleB_part3.FASTQ > sampleB.fastq
echo "sampleB.fastq created"

#Join Sample C files
##Concatemate the 3 parts for sample C into one complete file
cat sampleC_part1.FASTQ sampleC_part2.FASTQ sampleC_part3.FASTQ > sampleC.fastq
echo "sampleC.fastq created"

#Join Sample D files
##Concatemate the 3 parts for sample D into one complete file
cat sampleD_part1.FASTQ sampleD_part2.FASTQ sampleD_part3.FASTQ > sampleD.fastq
echo "sampleD.fastq created"

#Final Message
##All samples are now ready for FASTQC analysis
echo "All files concatenated"
