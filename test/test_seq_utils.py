from Bio.Seq import Seq
import unittest

from SeqUtils.seq_utils import read_fasta

class TestSeqUtils(unittest.TestCase):
    def test_read_fasta(self):
        """
        Test read_fasta function.
        """
        fasta_file = "test/test.fasta"
        sequences = read_fasta(fasta_file, output_type="dict")
        assert len(sequences) == 1
        assert sequences["seq_id1"] == Seq("AATCGGGAGAGATTCCCA")
        sequences = read_fasta(fasta_file, output_type="list")
        assert len(sequences) == 1
        assert sequences[0] == Seq("AATCGGGAGAGATTCCCA")
        sequences = read_fasta(fasta_file, output_type="pandas")
        assert len(sequences) == 1
        assert sequences["Sequence"][0] == Seq("AATCGGGAGAGATTCCCA")
        assert sequences["Description"][0] == "seq_id1 123 ACTG", sequences["Description"][0]
        assert sequences["ID"][0] == "seq_id1"
        
if __name__ == "__main__":
    unittest.main()
