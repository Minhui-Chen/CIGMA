#!/bin/bash
#SBATCH --job-name=AMM
#SBATCH --output=logs/AMM_%A_%a.out
#SBATCH --error=logs/AMM_%A_%a.err
#SBATCH --array=1-88
#SBATCH --time=10-00:00:00
#SBATCH --mem=64G
#SBATCH --partition=tier3q
#SBATCH --cpus-per-task=1


source activate amm

python AMM.py\
  --m 4\
  --iterator ${SLURM_ARRAY_TASK_ID}\
  --out ./\
  --pk_size 1 1 3 5 10 10 10 10\
  --set_names ./gene_set/set_name_norandom.txt
