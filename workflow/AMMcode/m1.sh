#!/bin/bash
#SBATCH --job-name=AMM
#SBATCH --output=logs/AMM_%A_%a.out
#SBATCH --error=logs/AMM_%A_%a.err
#SBATCH --array=1-22
#SBATCH --time=2-00:00:00
#SBATCH --mem=64G
#SBATCH --partition=tier3q
#SBATCH --cpus-per-task=1

source activate amm

python AMM.py \
	--m 1 \
	--iterator ${SLURM_ARRAY_TASK_ID} \
	--genes_ref ./cigma_gene_location/gnomad_gene_location_guide_cigma_chr \
	--snp_ref_bim /gpfs/data/ukb-share/dahl/minhuic/GxCTMM/data/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC. \
	--kn_k 50 \
	--snp_loc_col 3 \
	--genes_loc_col 1 \
	--out ./knmatrix_cigma/ \
	--snp_rsid_col 1 \
	--genes_ensg_col 0


