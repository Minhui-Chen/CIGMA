library(data.table)
data = read.csv('cis.jk.20250915.csv')
loc = fread('gene_location.txt')

cigma_gene = data$gene
loc = loc[which(loc$feature %in% cigma_gene),]
loc$mid = (loc$start + loc$end) / 2
for (i in 1:22) {
	a = loc[which(loc$chr == i),c('feature','mid')]
	outpath = paste0('./cigma_gene_location/gnomad_gene_location_guide_cigma_chr',i,'.txt')
	write.table(a, outpath, sep = ' ', row.names = FALSE, col.names = FALSE,
		quote = FALSE)
}

# remove MHC region genes from chr6
data = fread('gnomad_gene_location_guide_cigma_chr6.txt')
loc = fread('./gene_location.txt')
MHC_gene = loc[which(loc$chr == 6 & loc$start < 33448354 & loc$end > 28477797),]$'feature'

d = data[-which(data$V1 %in% MHC_gene),]
write.table(d, 'gnomad_gene_location_guide_cigma_chr6.txt', row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ' ')
