import unittest
import os
import subprocess

# run msa.R
msa_r_script = "../../Rscripts/msa.R"
inputpath="../test_data/test_seqs.csv"
outputpath="../test_data/test_seqs_msa.out.csv"
expected_file="../test_data/test_seqs_msa.expected.csv"
cmd = f"Rscript {msa_r_script} --input {inputpath} --output {outputpath} --group_index group --sequence sequence --seq_index seq_id --type dna"
subprocess.run(cmd.split())
with open(outputpath, 'r') as outfile, open(expected_file, 'r') as expected_file:
    assert outfile.read() == expected_file.read()