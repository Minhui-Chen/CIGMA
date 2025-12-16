
echo "./gene_set/shared_set" > shared.txt
echo "./gene_set/specific_set" > specific.txt
echo "./gene_set/DEG_set" > DEG.txt
echo "./gene_set/GCTA_set" > GCTA.txt


python AMM.py \
	--m 7 \
	--out ./ \
	--pk_size 1 1 3 5 10 10 10 10 \
	--control_name baseline_ld \
	--pk_out_name ./shared_alltraits.txt \
	--set_names ./shared.txt \
	--ss_list ./amm_gwas_sum.txt


python AMM.py \
	--m 7 \
	--out ./ \
	--pk_size 1 1 3 5 10 10 10 10 \
	--control_name baseline_ld \
	--pk_out_name ./specific_alltraits.txt \
	--set_names ./specific.txt \
	--ss_list ./amm_gwas_sum.txt


python AMM.py \
	--m 7 \
	--out ./ \
	--pk_size 1 1 3 5 10 10 10 10 \
	--control_name baseline_ld \
	--pk_out_name ./DEG_alltraits.txt \
	--set_names ./DEG.txt \
	--ss_list ./amm_gwas_sum.txt


python AMM.py \
	--m 7 \
	--out ./ \
	--pk_size 1 1 3 5 10 10 10 10 \
	--control_name baseline_ld \
	--pk_out_name ./GCTA_alltraits.txt \
	--set_names ./GCTA.txt \
	--ss_list ./amm_gwas_sum.txt
