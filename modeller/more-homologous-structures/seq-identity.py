#!/usr/bin/env python3.12

from Bio import SeqIO
from Bio import pairwise2

def read_fasta(file_path):
    """Read a FASTA file and return the sequence as a string."""
    with open(file_path, "r") as file:
        record = SeqIO.read(file, "fasta")
        return str(record.seq)

def calculate_sequence_identity(seq1, seq2):
    """Calculate the percentage of sequence identity between two sequences."""
    # Perform global pairwise alignment
    alignments = pairwise2.align.globalxx(seq1, seq2)
    # Get the best alignment
    best_alignment = alignments[0]
    
    # Extract alignment information
    aligned_seq1 = best_alignment[0]
    aligned_seq2 = best_alignment[1]
    
    # Calculate the number of matches
    matches = sum(a == b for a, b in zip(aligned_seq1, aligned_seq2) if a != '-' and b != '-')
    
    # Calculate the percentage of identity
    identity_percentage = (matches / len(seq1)) * 100
    
    return identity_percentage

# File paths to the FASTA files
file_path1 = "rcsb_pdb_7AAQ.fasta"
file_path2 = "glut9.fasta"

# Read sequences from FASTA files
seq1 = read_fasta(file_path1)
seq2 = read_fasta(file_path2)

# Calculate sequence identity
identity = calculate_sequence_identity(seq1, seq2)

# Print the result
print(f"Percentage sequence identity: {identity:.2f}%")
