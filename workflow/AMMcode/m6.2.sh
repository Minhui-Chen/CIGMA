#!/bin/bash
#SBATCH --job-name=AMM
#SBATCH --output=logs/AMM_%A_%a.out
#SBATCH --error=logs/AMM_%A_%a.err
#SBATCH --array=1-40
#SBATCH --time=10-00:00:00
#SBATCH --mem=64G
#SBATCH --partition=tier3q
#SBATCH --cpus-per-task=1


source activate amm

python AMM.py \
	--m 6 \
	--which_regression 2 \
	--iterator ${SLURM_ARRAY_TASK_ID} \
	--set_names ./gene_set/set_name_norandom.txt \
	--ss_list ./amm_gwas_sum.txt \
	--control_list ./baseline.txt \
	--weights ./data/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
	--freq_file ./data/ldsc/1000G_Phase3_frq/1000G.EUR.QC. \
	--ldsc_path ./ldsc/ldsc.py \
	--out ./

