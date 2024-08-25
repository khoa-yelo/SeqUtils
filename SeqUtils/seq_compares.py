"""
Module to compare sequence distances and similarities.
Khoa Hoang
08/06/2024
"""

import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
from polyleven import levenshtein
from seq_utils import read_fasta
import fire

def hamming_distance(seq1, seq2):
    """
    Calculate the Hamming distance between two sequences of equal length.

    Parameters:
    - seq1: The first input sequence.
    - seq2: The second input sequence.

    Returns:
    - The Hamming distance between the two sequences.
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length.")
    
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def jaccard_similarity(set1, set2):
    """
    Calculate the Jaccard similarity between two sets.

    Parameters:
    - set1: The first input set.
    - set2: The second input set.

    Returns:
    - The Jaccard similarity between the two sets.
    """
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    
    return intersection / union

def levenshtein_distance(seq1, seq2):
    """
    Calculate the Levenshtein distance between two sequences.

    Parameters:
    - seq1: The first input sequence.
    - seq2: The second input sequence.

    Returns:
    - The Levenshtein distance between the two sequences.
    """
    # n, m = len(seq1), len(seq2)
    # dp = np.zeros((n + 1, m + 1))
    
    # for i in range(n + 1):
    #     dp[i][0] = i
    # for j in range(m + 1):
    #     dp[0][j] = j
    
    # for i in range(1, n + 1):
    #     for j in range(1, m + 1):
    #         cost = 0 if seq1[i - 1] == seq2[j - 1] else 1
    #         dp[i][j] = min(dp[i - 1][j] + 1, dp[i][j - 1] + 1, dp[i - 1][j - 1] + cost)

    dist = levenshtein(seq1, seq2)
    
    return dist


def compare_sequences(seq1, seq2, method="hamming"):
    """
    Compare two sequences using a specified method.

    Parameters:
    - seq1: The first input sequence.
    - seq2: The second input sequence.
    - method: The method to use for comparison (default: "hamming").

    Returns:
    - The distance or similarity between the two sequences.
    """
    if method == "hamming":
        return hamming_distance(seq1, seq2)
    elif method == "jaccard":
        return jaccard_similarity(set(seq1), set(seq2))
    elif method == "levenshtein":
        return levenshtein_distance(seq1, seq2)
    else:
        raise ValueError("Invalid method. Choose from 'hamming', 'jaccard', or 'levenshtein'.")


def compare_sequences_fasta(fasta_file, method="hamming", output_file=None):
    """
    Compare sequences in a fasta file using a specified method.

    Parameters:
    - fasta_file: The input fasta file containing sequences.
    - method: The method to use for comparison (default: "hamming").
    - output_file: The output file to write the results to.
    """
    results = []
    SLURM_ARRAY_TASK_ID = os.getenv("SLURM_ARRAY_TASK_ID")
    SLURM_ARRAY_TASK_COUNT = os.getenv("SLURM_ARRAY_TASK_COUNT")
    # use tqdm to show progress bar
    records = read_fasta(fasta_file, output_type="dict")
    print("Total sequences in fasta file: ", len(list(records)))

    if SLURM_ARRAY_TASK_ID is None:
        print("SLURM_ARRAY_TASK_ID is not set. Running on all sequences.")
        for key1, seq1 in records.items():
            for key2, seq2 in records.items():
                result = compare_sequences(str(seq1), str(seq2), method)
                results.append((key1, key2, result))
    else:
        i = 0
        for key1, seq1 in records.items():
            if i % int(SLURM_ARRAY_TASK_COUNT) == int(SLURM_ARRAY_TASK_ID):
                for key2, seq2 in records.items():
                    result = compare_sequences(str(seq1), str(seq2), method)
                    results.append((key1, key2, result))
            i += 1

    df = pd.DataFrame(results, columns=["Sequence1", "Sequence2", method.capitalize()])
    df = df.pivot(index="Sequence1", columns="Sequence2", values=method.capitalize())
    if output_file:
        output_file = output_file.replace(".csv", f"_{SLURM_ARRAY_TASK_ID}.csv") if SLURM_ARRAY_TASK_ID else output_file
        df.to_csv(output_file, index=False)
    else:
        print(df)


if __name__ == "__main__":
    # seq1 = "ATTCCATCATATTCA"
    # seq2 = "ATTCC_ATCATATTCC"
    # method = "levenshtein"
    # result = compare_sequences(seq1, seq2, method)
    # print(f"{method.capitalize()} distance between the two sequences: {result}")
    # fasta_file = "/oak/stanford/groups/horence/khoa/scratch/repos/SeqUtils/test/test_data/test.fasta"
    # output_file = "/oak/stanford/groups/horence/khoa/scratch/repos/SeqUtils/test/test_data/test_levdist.csv"
    # compare_sequences_fasta(fasta_file, method, output_file)
    fire.Fire()