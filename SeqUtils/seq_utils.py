import Bio 
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import os
from os.path import join
import sys

SCRATCH = os.getenv("SCRATCH")
REPO = join(SCRATCH, "repos", "SeqUtils")

def read_fasta(fasta_file, output_type="dict"):
    """
    Read a fasta file and return a dictionary with the sequence id as key and the sequence as value.
    """
    if output_type == "list":
        sequences = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(record.seq)
        return sequences
    elif output_type == "dict":
        sequences = {}
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences[record.id] = record.seq
        return sequences
    elif output_type == "pandas":
        sequences = []
        description = []
        ids = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(record.seq)
            description.append(record.description)
            ids.append(record.id)
        return pd.DataFrame({"ID": ids, "Description": description, "Sequence": sequences})


def write_fasta(df, output_file, id_col, description_col, sequence_col):
    """
    Write a pandas dataframe to a fasta file.
    """
    with open(output_file, "w") as f:
        for index, row in df.iterrows():
            f.write(">" + str(row[id_col]) + " " + str(row[description_col]) + "\n")
            f.write(str(row[sequence_col]) + "\n")


def split_fasta(fasta_file, output_dir):
    """
    Split a fasta file into multiple files.
    """
    for record in SeqIO.parse(fasta_file, "fasta"):
        output_file = os.path.join(output_dir, record.id + ".fasta")
        with open(output_file, "w") as f:
            f.write(">" + record.description + "\n")
            f.write(str(record.seq) + "\n")