#importing the module I need 
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
import os
import csv
import sys

# Fasta file path.
fasta_infile = "protein_sequences.fasta"

# Iterate over all sequences within the fasta file.
for record in SeqIO.parse(fasta_infile, "fasta"):
    seq_id = record.id
    seq_description = record.description
    
    #the amino acids left over from the cut site
    leftoverres = ("SN")

    # sequnece of the construct with all noncritical residues
    protein_seq_record = record.seq

    #the amount of non-critical residues on the c terminal
    input_CterminalNC = 55

    # the amount of non-critical residues on the N terminal
    input_NterminalNC = 57

    CterminalNC = int(input_CterminalNC)
    NterminalNC = int(input_NterminalNC)

    protein_seq = str(protein_seq_record)

    #csv_outfile = os.path.join(output_dir, "test.csv")
    basename = os.path.basename(fasta_infile)
    filename = os.path.splitext(basename)[0]
    
    csv_outfile = "_".join([seq_id, filename]) + ".tsv"
    csv_output_file = open(csv_outfile, 'w+')
    csv_writer = csv.writer(csv_output_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_NONE)

    # Will be for header. Can make it work by getting the numbers before hand.
    #csv_writer.writerow([])

    # Get the length of the protein sequence so we can remove amino acids from the C terminus end.
    protein_length = len(protein_seq)
    
    # Starting at the C terminal region of the protein at the length of the protein.
    for CterminalNC_value in range(protein_length, (protein_length - CterminalNC) - 1, -1):
        
        c1list = []
        
        # Protein sequence starting with full length.
        seq = ""
        for input_NterminalNC_value in range(0, input_NterminalNC, 1):
            
                        
            seq = protein_seq[input_NterminalNC_value:CterminalNC_value]
            
            analysed_seq = ProteinAnalysis(leftoverres + seq)
            analysed_seq.instability_index()
            print(", ".join([str(input_NterminalNC_value),str(CterminalNC_value), str(round(analysed_seq.instability_index(), 2))]))
            print(seq)
            c1list.append(round(analysed_seq.instability_index(), 2))

        csv_writer.writerow(c1list)


