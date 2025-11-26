setwd('./target/')

Ds = c('CAD','SCZ','Celiac','Crohns','SLE','MS','PBC','RA','UC')
cigma = read.csv('./data/cis.jk.20250915.csv')
gene_info = fread('all_humangene_ncbi.txt')
cigma_noMHC = fread('../../cigma_gene_noMHC10178.txt')
cigma = cigma[which(cigma$gene %in% cigma_noMHC$x),]

col = 'free.v'  # or
col = 'free.hom_g2'

res_mean = matrix(NA, nrow = length(Ds), ncol = 6)

for (i in 1:length(Ds)) {
	trait = Ds[i]
	data = fread(paste0(trait,'.tsv'))
	if (sum(is.na(data$chembl)) != 0) { print('ERROR') }
	targets = data$symbol
	row1 = which(gene_info$Symbol %in% targets)
	left_targets = targets[-which(targets %in% gene_info$Symbol)]
	row2 = which(sapply(gene_info$Synonyms, function(x) {
    	any(strsplit(x, "\\|")[[1]] %in% left_targets) }))
	row = union(row1,row2)
	gene_info_subset <- gene_info[row, ]
	gene_info_subset = gene_info_subset[which(!is.na(gene_info_subset$Ensembl_id)),]
	target = gene_info_subset$Ensembl_id

	cigma_target = cigma[which(cigma$gene %in% target),]
	cigma_nontarget = cigma[which(!cigma$gene %in% target),]

	res_mean[i,1] = Ds[i]
	t = t.test(cigma_target[,col],cigma_nontarget[,col])

	res_mean[i,2] = t$p.value
	res_mean[i,3] = mean(cigma_target[,col])
	res_mean[i,4] = sqrt(var(cigma_target[,col]) / length(cigma_target[,col]))
	res_mean[i,5] = mean(cigma_nontarget[,col])
	res_mean[i,6] = sqrt(var(cigma_nontarget[,col]) / length(cigma_nontarget[,col]))
}

res_mean = as.data.frame(res_mean)
colnames(res_mean) = c('disease','pval','mean_target','se_mean_target','mean_nontarget','se_mean_nontarget')
res_mean

