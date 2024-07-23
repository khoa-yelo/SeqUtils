"""Module containing functions for counting biological sequences."
Khoa Hoang
07/22/2024
"""

import os
from Bio.Seq import Seq
from Bio import SeqIO
import pandas as pd


def get_kmers(sequence, k, i):
    """
    Get k-mers of length k from the sequence starting at every ith index.

    Parameters:
    - sequence: The input string or sequence from which to extract k-mers.
    - k: The length of each k-mer.
    - i: The step size to obtain k-mers at every ith index.

    Returns:
    - A set of k-mers.
    """
    kmers = []
    
    # Loop through the sequence, starting at each ith index
    for start in range(0, len(sequence) - k + 1, i):
        kmer = sequence[start:start + k]
        kmers.append(kmer)
    
    return list(set(kmers))

def get_kmers_fasta(fasta_file, k, i, output_file):
    """
    Get k-mers of length k from sequences in a fasta file starting at every ith index.

    Parameters:
    - fasta_file: The input fasta file containing sequences.
    - k: The length of each k-mer.
    - i: The step size to obtain k-mers at every ith index.
    """
    all_kmers = []    
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        kmers = get_kmers(seq, k, i)
        all_kmers.extend(kmers)
    all_kmers = set(all_kmers)
    with open(output_file, "w") as f:
        for kmer in all_kmers:
            f.write(kmer + "\n")

if __name__ == "__main__":
    sequence = "ATTCCATCTATTCA"
    k = 5
    i = 3
    kmers = get_kmers(sequence, k, i)
    print(kmers)