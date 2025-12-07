from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# input and output files
input_file = "combined_all_seq.fasta"
output_file = "AAtranslated_seqs.fasta"

# read the DNA sequences
print(f"Reading: {input_file}")
sequences = list(SeqIO.parse(input_file, "fasta"))
print(f"Loaded {len(sequences)} sequences")
print()

# list to store translates AA sequences
aa_sequences = []
# translate each sequence
for record in sequences:
    dna_seq = record.seq
    best_orf = ""
    best_frame = None
    for frame in range(3):
        # Table 2 for vertebrate mitochondrial
        protein = dna_seq[frame:].translate(table=2, to_stop=False)
        fragments = str(protein).split("*")
        longest_in_frame = max(fragments, key=len)
        if len(longest_in_frame) > len(best_orf):
            best_orf = longest_in_frame
            best_frame = frame
    best_protein = Seq(best_orf)
    aa_record = SeqRecord(
        best_protein,
        id=record.id,
        name=record.name,
        description=f"{record.description} [frame={best_frame}]",
        annotations={
            "type": "longest ORF",
            "table": 2,
            "frame": best_frame
        }
    )

    aa_sequences.append(aa_record)

# Print progress
    print(f"{record.id:25s} {len(dna_seq):5d} bp → {len(best_protein):4d} aa (frame {best_frame})")
print()

# write translated sequences to FASTA file
SeqIO.write(aa_sequences, output_file, "fasta")
