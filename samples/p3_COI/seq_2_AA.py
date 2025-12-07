from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#File Set Up
##Input and Output files
input_file = r"c:/Users/gt25963/meat_species/samples/p3_COI/FINAL_comb_all_seq.fasta"
output_file = "FINAL_AA_translated_seqs.fasta"

#Load in the DNA sequences from the FASTA file
print(f"Reading: {input_file}")
sequences = list(SeqIO.parse(input_file, "fasta"))
print(f"Loaded {len(sequences)} sequences")
print()

#List to store oputput amino-acid SeqRecord objects
aa_sequences = []
#Translate each sequence
for record in sequences:
    dna_seq = record.seq #Nucleotide sequence from the FASTA record
    best_orf = "" #Will store the longest ORF found across all frames
    best_frame = None #Will store which frame (0,1,2) contains the longest ORF
    #Translate sequence in all three forward reading frames
    for frame in range(3):
        #Translate starting at this frame using vertebrate mitochondiral code (table=2). Keeping stop codons (to_stop=Flase) as to identify ORFs
        protein = dna_seq[frame:].translate(table=2, to_stop=False) 
        #Split translation into fragments separated by stop codons ("*"). Each fragment is a potential ORF
        fragments = str(protein).split("*")
        #Find the longest continuous amino-acid stretch (ORF) in this frame
        longest_in_frame = max(fragments, key=len)
        #If this ORF is longer than the best one seen so far, store it.
        if len(longest_in_frame) > len(best_orf):
            best_orf = longest_in_frame
            best_frame = frame
    #Convert the longest ORF strong back into Seq object
    best_protein = Seq(best_orf)
    #Create a new SeqRecord containing the longest ORF for this sequence
    aa_record = SeqRecord(
        best_protein,
        id=record.id, #Keep the same FASTA id
        name=record.name, #Keep the original name
        description=f"{record.description} [frame={best_frame}]", #Add frame information into description
        annotations={
            "type": "longest ORF", #Mark this as the longest ORF
            "table": 2, #Used the vertebrate mitochondrial code (table 2)
            "frame": best_frame #Which reading frame the ORF came from
        }
    )

    aa_sequences.append(aa_record) #Add to output list

#Print the progress for each sequence 
    print(f"{record.id:25s} {len(dna_seq):5d} bp -> {len(best_protein):4d} aa (frame {best_frame})")
print()

#Write all translated longest-ORF sequences to a FASTA file
SeqIO.write(aa_sequences, output_file, "fasta")
