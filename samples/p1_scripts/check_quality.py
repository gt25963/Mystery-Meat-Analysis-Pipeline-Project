from Bio import SeqIO

samples = ['sampleA.fastq', 'sampleB.fastq', 'sampleC.fastq', 'sampleD.fastq']

for sample_file in samples:
    print(f"\n=== {sample_file} ===")
    try:
        for record in SeqIO.parse(sample_file, "fastq"):
            avg_quality = sum(
                record.letter_annotations['phred_quality']) / len(record)
            print(f"Average Phred score: {avg_quality:.2f}")
    except ValueError as e:
        print(f"ERROR in {sample_file}: {e}")
        print("This file may have formatting issues - check for spaces in sequence lines")
