import os
from os.path import join
from Bio.Seq import Seq
import unittest
import pandas as pd
from SeqUtils.seq_utils import read_fasta, write_fasta, split_fasta, REPO
from SeqUtils.seq_mod import back_translate
import sys


class TestSeqUtils(unittest.TestCase):
    def test_read_fasta(self):
        """
        Test read_fasta function.
        """
        fasta_file = "test/test_data/test.fasta"
        sequences = read_fasta(fasta_file, output_type="dict")
        assert len(sequences) == 2
        assert sequences["seq_id1"] == Seq("AATCGGGAGAGATTCCCA")
        sequences = read_fasta(fasta_file, output_type="list")
        assert len(sequences) == 2
        assert sequences[0] == Seq("AATCGGGAGAGATTCCCA")
        sequences = read_fasta(fasta_file, output_type="pandas")
        assert len(sequences) == 2
        assert sequences["Sequence"][0] == Seq("AATCGGGAGAGATTCCCA")
        assert sequences["Description"][0] == "seq_id1 123 ACTG", sequences["Description"][0]
        assert sequences["ID"][0] == "seq_id1"

    def test_write_fasta(self):
        """
        Test write_fasta function.
        """
        csv_file = "test/test_data/test_seqs.csv"
        df = pd.read_csv(csv_file)
        output_file = "test/test_data/test_write_fasta.out.fasta"
        expected_file = "test/test_data/test_write_fasta.expected.fasta"
        write_fasta(df, output_file, id_col="seq_id", description_col="group", sequence_col="sequence")
        with open(output_file, 'r') as outfile, open(expected_file, 'r') as expected_file:
            assert outfile.read() == expected_file.read()
    
    def test_back_translate(self):
        """
        Test back_translate function.
        """
        protein_seq = "MAG"
        dna_seq = back_translate(protein_seq)
        assert dna_seq == "ATGGCTGGT", dna_seq

    def test_split_fasta(self):
        """
        Test split_fasta function.
        """
        fasta_file = join(REPO, "test/test_data/test.fasta")
        output_dir = join(REPO, "test/test_data/split_fasta")
        os.makedirs(output_dir, exist_ok=True)
        split_fasta(fasta_file, output_dir)
        assert os.path.exists(join(output_dir, "seq_id1.fasta"))
        assert os.path.exists(join(output_dir, "seq_id2.fasta"))
        
if __name__ == "__main__":
    unittest.main()
