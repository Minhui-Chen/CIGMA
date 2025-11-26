for trait in Height CAD SCZ UC RA PBC MS Crohns Celiac Lupus
do
	python AMM.py \
	--m 8 \
	--out ./ \
	--set_names ./gene_set/set_name_norandom.txt \
	--control_name baseline_ld \
	--control_name_m /gpfs/data/ukb-share/dahl/minhuic/GxCTMM/data/ldsc/1000G_Phase3_baselineLD_v2.2/baselineLD \
	--ss_list ./$trait.txt \
	--n_genes_total 10178

	awk 'FNR==1 && NR!=1 { next } { print }' ./gene_set/shared_set.baseline_ld.enrichment ./gene_set/specific_set.baseline_ld.enrichment ./gene_set/DEG_set.baseline_ld.enrichment ./gene_set/GCTA_set.baseline_ld.enrichment > enrichment.$trait.txt
done

awk 'FNR==1 && NR!=1 { next } { print }' enrichment.*.txt > enrichment.merged.txt
