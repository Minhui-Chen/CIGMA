#!/bin/bash
#SBATCH --job-name=AMM
#SBATCH --output=logs/AMM_%A_%a.out
#SBATCH --error=logs/AMM_%A_%a.err
#SBATCH --array=1-16
#SBATCH --time=10-00:00:00
#SBATCH --mem=64G
#SBATCH --partition=tier3q
#SBATCH --cpus-per-task=1

source activate amm

python AMM.py \
  --m 2 \
  --iterator ${SLURM_ARRAY_TASK_ID} \
  --set_names ./gene_set/set_name.txt \
  --kn_in_dir ./knmatrix_cigma/ \
  --out ./ \
  --kn_k 50 \
  --kn_k_out 50