# shared
for i in {1..22}
do
  python AMM.py \
	--m 5 \
	--iterator $i \
	--set_names ./shared.txt \
	--pk_vector 0.146152459 0.411924851 -0.213547788 -0.101820128 -0.028405088 -0.034737604 0.141205095 0.081104266 \
	--pk_size 1 1 3 5 10 10 10 10 \
	--out ./
done

# specific
for i in {1..22}
do
  python AMM.py \
	--m 5 \
	--iterator $i \
	--set_names ./specific.txt \
	--pk_vector 0.065569724 0.080138792 0.050650402 -0.0007467 0.016788946 0.036747877 0.011036483 0.006034072 \
	--pk_size 1 1 3 5 10 10 10 10 \
	--out ./
done


# DEG
for i in {1..22}
do
  python AMM.py \
	--m 5 \
	--iterator $i \
	--set_names ./DEG.txt \
	--pk_vector 0.254919955 0.048141104 0.199694491 -0.03282177 0.033974925 0.000816415 0.008702587 -0.017297495 \
	--pk_size 1 1 3 5 10 10 10 10 \
	--out ./
done

# GCTA
for i in {1..22}
do
  python AMM.py \
	--m 5 \
	--iterator $i \
	--set_names ./GCTA.txt \
	--pk_vector -0.042459257 0.009632092 -0.538958297 -0.052682327 0.107734203 0.02289006 0.175115202 -0.014428096 \
	--pk_size 1 1 3 5 10 10 10 10 \
	--out ./
done

