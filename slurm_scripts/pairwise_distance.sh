#!/bin/bash
#SBATCH --job-name=pairwise_distance
#SBATCH --output=pairwise_distance.%j.out
#SBATCH --error=pairwise_distance.%j.err
#SBATCH --time=02:10:00
#SBATCH -p normal,horence
#SBATCH -c 1
#SBATCH --mem=8GB
#SBATCH --array=0-9


python3 /oak/stanford/groups/horence/khoa/scratch/repos/SeqUtils/SeqUtils/seq_compares.py compare_sequences_fasta \
--fasta_file /oak/stanford/groups/horence/khoa/scratch/repos/RestrictionEnzyme/data/all_enzymes_back_translated.fasta \
--method levenshtein \
--output_file /oak/stanford/groups/horence/khoa/scratch/repos/RestrictionEnzyme/data/all_enzymes_back_translated_lev_dist.csv
