"""Module containing functions for modifying biological sequences.
Khoa Hoang
07/22/2024
"""

import os
from Bio.Seq import Seq
import pandas as pd
from Bio import SeqIO

CODON_TABLE = {
    'A': 'GCT', 'C': 'TGT', 'D': 'GAT', 'E': 'GAA', 'F': 'TTT',
    'G': 'GGT', 'H': 'CAT', 'I': 'ATT', 'K': 'AAA', 'L': 'CTT',
    'M': 'ATG', 'N': 'AAT', 'P': 'CCT', 'Q': 'CAA', 'R': 'CGT',
    'S': 'TCT', 'T': 'ACT', 'V': 'GTT', 'W': 'TGG', 'Y': 'TAT',
}

def back_translate(protein_seq):
    """back translate from DNA to protein"""
    dna_seq = []
    for aa in protein_seq:
        codon = CODON_TABLE.get(aa)
        if codon is None:
            raise ValueError(f"No codon found for amino acid: {aa}")
        dna_seq.append(codon)
    return ''.join(dna_seq)


def back_translate_fasta(fasta_file, output_file):
    """
    Back translate a fasta file.
    """
    with open(output_file, "w") as f:
        for record in SeqIO.parse(fasta_file, "fasta"):
            dna_seq = back_translate(str(record.seq))
            f.write(">" + record.description + "\n")
            f.write(dna_seq + "\n")


def reverse_complement_fasta(fasta_file, output_file):
    """
    Reverse complement a fasta file.
    """
    with open(output_file, "w") as f:
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq = record.seq.reverse_complement()
            f.write(">" + record.description + "\n")
            f.write(str(seq) + "\n")


def reverse_fasta(fasta_file, output_file):
    """
    Reverse a fasta file.
    """
    with open(output_file, "w") as f:
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq = record.seq[::-1]
            f.write(">" + record.description + "\n")
            f.write(str(seq) + "\n")