# Mystery Meat Analysis Pipeline

A bioinformatics pipeline for identifying unknown meat samples using DNA barcoding and phylogenetic analysis of the COI (Cytochrome c Oxidase subunit I) gene. This project demonstrates the complete workflow from raw sequencing data to phylogenetic tree construction, using quality control, BLAST searches, sequence translation, and alignment techniques.

## Repository Structure

| Directory/File | Description |
|----------------|-------------|
| `samples/p1_scripts/` | Phase 1: Sequence concatenation and quality control |
| `samples/p1_scripts/concat_FASTQC.sh` | Main script to concatenate all four samples |
| `samples/p1_scripts/concatA_fastq.sh` | Concatenate sample A parts |
| `samples/p1_scripts/concatB_fastq.sh` | Concatenate sample B parts |
| `samples/p1_scripts/concatC_fastq.sh` | Concatenate sample C parts |
| `samples/p1_scripts/concatD_fastq.sh` | Concatenate sample D parts |
| `samples/p1_scripts/check_quality.py` | Quality assessment script |
| `samples/p2_COI/` | Phase 2: Reference COI sequences from marine species |
| `samples/p2_COI/concat_ref_seqs.sh` | Combine all reference sequences |
| `samples/p3_COI/` | Phase 3: Translation and alignment |
| `samples/p3_COI/seq_2_AA.py` | Translate DNA to amino acids (identifies longest ORF) |
| `requirements.txt` | Python package dependencies |
| `.gitignore` | Git ignore rules for large data files |
| `README.md` | Project documentation |

## Features

* Concatenates split FASTQ files from sequencing runs
* Quality control of raw sequencing data
* BLAST searches against COI reference database
* Identifies longest open reading frame (ORF) across three reading frames
* Translates DNA sequences to proteins using vertebrate mitochondrial genetic code (table 2)
* Multiple sequence alignment of translated sequences
* Phylogenetic analysis of mystery meat samples against reference species
* Reference database includes marine mammals (whales, dolphins), sea turtles, and tuna species

## Installation

Install required Python packages:

```bash
pip install biopython
pip install pandas
pip install matplotlib
pip install numpy
```

Or install all dependencies at once:

```bash
pip install -r requirements.txt
```

Additional software required:
* **FastQC** - for quality control
* **BLAST+** - for sequence searching
* **Clustal Omega** - for multiple sequence alignment

## Usage

### Phase 1: Concatenate Sequences and Quality Control

Navigate to the scripts directory and concatenate FASTQ files:

```bash
cd samples/p1_scripts/
bash concat_FASTQC.sh
```

This combines split FASTQ files (part1, part2, part3) into complete sequences for each sample.

### Phase 2: Prepare Reference Database

Reference sequences are stored in `samples/p2_COI/` and include:
* Marine mammals: *Balaena mysticetus*, *Eubalaena glacialis*, *Delphinapterus leucas*, *Lagenorhynchus albirostris*
* Sea turtles: *Dermochelys coriacea*, *Caretta caretta*, *Eretmochelys imbricata*, *Natator depressus*
* Tuna: *Thunnus thynnus*, *Thunnus albacares*, *Thunnus maccoyii*, *Sarda sarda*

Combine all reference sequences:

```bash
cd ../p2_COI/
bash concat_ref_seqs.sh
```

### Phase 3: Translate and Align

Translate DNA sequences to amino acids using the longest ORF method:

```bash
cd ../p3_COI/
python seq_2_AA.py
```

This script:
* Translates sequences in all three reading frames
* Identifies the longest continuous amino acid sequence (ORF)
* Uses vertebrate mitochondrial genetic code (table 2)
* Outputs protein sequences ready for alignment

Align translated sequences using Clustal Omega:

```bash
clustalo -i FINAL_AA_translated_seqs.fasta -o FINAL_Alignment.fa --auto
```

## Output Files

* `sampleA.fastq`, `sampleB.fastq`, `sampleC.fastq`, `sampleD.fastq` - Concatenated sequencing reads
* `FINAL_all_ref_sequences.fasta` - Combined reference database
* `FINAL_AA_translated_seqs.fasta` - Translated protein sequences
* `FINAL_Alignment.fa` - Multiple sequence alignment
* `FINAL_alignment_image.png` - Alignment visualization

## Methodology

This pipeline uses DNA barcoding with the COI gene, a standard marker for animal species identification. The gene is highly conserved within species but variable between species, making it ideal for taxonomic classification.

The translation workflow identifies open reading frames by:
1. Translating DNA in all three forward reading frames
2. Splitting translations on stop codons to identify ORFs
3. Selecting the longest continuous amino acid sequence
4. Using vertebrate mitochondrial genetic code (table 2) appropriate for COI genes

## Author

**Martina**  
MSc Bioinformatics, University of Bristol  
CEO, SentiDx

## Educational Context

This pipeline was developed as coursework for MSc Bioinformatics at the University of Bristol, demonstrating skills in bash scripting, Python/Biopython programming, sequence analysis workflows, and phylogenetic methods.
