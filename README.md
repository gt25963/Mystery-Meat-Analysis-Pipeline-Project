# Mystery Meat Analysis Pipeline

This project is a bioinformatics pipeline for identifying unknown meat samples using DNA barcoding and phylogenetic analysis of the COI (Cytochrome c Oxidase subunit I) gene. This project demonstrates the complete workflow from raw sequencing data to phylogenetic tree construction, using quality control, BLAST searches, sequence translation, and alignment techniques.

## Repository Structure

| Directory/File | Description |
|----------------|-------------|
| `samples/p1_scripts/` | Part 1: Data Reformatting and Quality Control Scripts |
| `samples/p1_scripts/concat_FASTQC.sh` | Script to concatenate the split fastq files of the four samples for FASTQC analysis |
| `samples/p1_scripts/concatA_fastq.sh` | Concatenate parts 1,2,3 of sample A |
| `samples/p1_scripts/concatB_fastq.sh` | Concatenate parts 1,2,3 of sample B |
| `samples/p1_scripts/concatC_fastq.sh` | Concatenate parts 1,2,3 of sample C |
| `samples/p1_scripts/concatD_fastq.sh` | Concatenate parts 1,2,3 of sample D |
| `samples/p2_COI/` | Part 2: Reference COI Sequences from NCBI BLAST Searches |
| `samples/p2_COI/concat_ref_seqs.sh` | Combine all reference sequences |
| `samples/p3_COI/` | Part 3: Sequence Translation and Multiple Sequence Alignment |
| `samples/p3_COI/seq_2_AA.py` | Translate DNA to Amino Acids (identifies longest ORF) |
| `README.md` | Project documentation |

## Features

* Concatenates split FASTQ files
* Quality control of raw sequencing data
* BLAST searches against NCBI nucelotide database (COI regions)
* Identifies longest open reading frame (ORF) across three reading frames
* Translates DNA sequences to proteins using vertebrate mitochondrial genetic code (table 2)
* Multiple sequence alignment of translated sequences
* Phylogenetic analysis of mystery meat samples against reference species
* Reference database includes marine mammals (whales, dolphins), sea turtles, and tuna species

## Installation

Install required Python packages:

```bash
pip install biopython
```

Additional software required:
* **FastQC** - for quality control
* **BLAST+** - for sequence searching
* **Clustal Omega** - for multiple sequence alignment
* **Jalview** - for alignment visualization (optional)

## Usage

### Part 1: Concatenate Sequences and Quality Control

For FASTQC analysis:

```bash
cd samples/p1_scripts/
bash concat_FASTQC.sh
```

For BLAST searches:

```bash
cd samples/p1_scripts/
bash concatA_fastq.sh
```

Combines split FASTQ files (part1, part2, part3) into complete sequences for each sample.

### Part 2: Prepare Reference Database

Reference sequences are stored in `samples/p2_COI/` and include:
* Marine mammals: *Balaena mysticetus*, *Eubalaena glacialis*, *Eubalaena japonica*, *Eschrichtius robustus*, *Delphinapterus leucas*, *Lagenorhynchus albirostris*, *Cephalorhynchus eutropia*, *Cephalorhynchus heavisidii*
* Sea turtles: *Dermochelys coriacea*, *Caretta caretta*, *Eretmochelys imbricata*, *Natator depressus*
* Tuna: *Thunnus thynnus*, *Thunnus albacares*, *Thunnus maccoyii*, *Sarda sarda*

This pipeline uses DNA barcoding with the COI gene, a standard for animal species identification. 

Combine all reference sequences:

```bash
cd ../p2_COI/
bash concat_ref_seqs.sh
```

### Part 3: Translate DNA to Amino Acids and Align the Sequences

Translate DNA sequences to amino acids using the longest ORF method:

```bash
cd ../p3_COI/
python seq_2_AA.py
```

Align translated sequences using Clustal Omega on the EBI web interface.

## Output Files

* `sampleA.fastq`, `sampleB.fastq`, `sampleC.fastq`, `sampleD.fastq` - Concatenated sample sequences
* `FINAL_all_ref_sequences.fasta` - Combined reference sequences
* `FINAL_AA_translated_seqs.fasta` - Translated protein sequences
* `FINAL_Alignment.fa` - Multiple sequence alignment
* `FINAL_alignment_image.png` - Alignment visualization
