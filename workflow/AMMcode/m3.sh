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
  --m 3\
  --iterator ${SLURM_ARRAY_TASK_ID} \
  --set_names ./gene_set/set_name_norandom.txt \
  --ldsc_path ./ldsc/ldsc.py \
  --lds_ref_binary ./data/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC. \
  --lds_snps_out ./snps_1000G_hapmap3/snps.\
  --out ./
