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
    left_over_res = ("SN")

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

    # Get the length of the protein sequence so we can remove amino acids from the C terminus end.
    protein_length = len(protein_seq)
    
    # Get the amino acid position for the N-terminal for the header of the output file.
    header_list = []
    for NterminalNC_index in range(0, NterminalNC, 1):
        #print((NterminalNC_index + 1))
        NterminalNC_position = (NterminalNC_index + 1)
        header_list.append(str(NterminalNC_position))
        
    # Write the header row.
    csv_writer.writerow([""] + header_list)
    
    # Delete header list.
    del header_list
    
    # Protein sequence starting with full length.
    seq = protein_seq
    
    # Starting at the C terminal region of the protein at the length of the protein.
    for CterminalNC_index in range(protein_length, (protein_length - CterminalNC), -1):
        
        instability_index_list = []
        
        CterminalNC_position = CterminalNC_index
        for NterminalNC_index in range(0, NterminalNC, 1):
            
            seq = protein_seq[NterminalNC_index:CterminalNC_index]
            analysed_seq = ProteinAnalysis(left_over_res + seq)
            analysed_seq.instability_index()
            
            print(", ".join([str(NterminalNC_index),str(CterminalNC_index), str(round(analysed_seq.instability_index(), 2))]))
            print(seq)
            
            instability_index_list.append(round(analysed_seq.instability_index(), 2))
            


        
        csv_writer.writerow([str(CterminalNC_position)] + instability_index_list)


