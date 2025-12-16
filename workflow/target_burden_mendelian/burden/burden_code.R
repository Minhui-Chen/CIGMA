setwd('./burden/')
library(data.table)

setwd('./mis/') # plof or syn
Ds = c('height','CAD','SCZ','Celiac','Crohns','SLE','MS','RA','UC')
cigma = read.csv('../cis.jk.csv')
cigma_noMHC = fread('../cigma_gene_noMHC10178.txt')
cigma = cigma[which(cigma$gene %in% cigma_noMHC$x),]


cols = c("free.v", "free.hom_g2")
sel_test = c("P.Value.SKATO", "P.Value.Burden", "P.Value.SKAT")


for (t in 1:length(sel_test)) {
	test = sel_test[t]
	
	res = matrix(NA, nrow = length(Ds)*length(cols), ncol = 7)

	for (c in 1:length(cols)) {
		col = cols[c]
		for (i in 1:length(Ds)) {
			trait = Ds[i]
			data = read.csv(paste0(trait,'.csv'))
			colnames(data)[2] = 'gene'
			d = merge(data, cigma, by='gene')

			thr = 1e-3
			d_burden = d[which(d[,test] < thr),]
			d_nonburden = d[which(!d[,test] < thr),]

			t = t.test(d_burden[,col], d_nonburden[,col])
			row_id = (c-1)*length(Ds) + i
			res[row_id , 1] = trait
			res[row_id , 2] = t$p.value
			res[row_id , 3] = mean(d_burden[,col])
			res[row_id , 4] = sqrt(var(d_burden[,col]) / length(d_burden[,col]))
			res[row_id , 5] = mean(d_nonburden[,col])
			res[row_id , 6] = sqrt(var(d_nonburden[,col]) / length(d_nonburden[,col]))
			res[row_id , 7] = col

		}
	}

	res = as.data.frame(res)
	colnames(res) = c('disease','t-test-P','mean_burden','se_burden','mean_nonburden','se_nonburden','variable')
	outpath = paste0(test,'_noMHC','.csv')
	write.csv(res, outpath, row.names = FALSE, quote = TRUE)
}





