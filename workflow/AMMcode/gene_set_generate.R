library(data.table)

res = read.csv('./cis.jk.csv')
loc = fread('./gene_location.txt')
sum(res$gene %in% loc$feature)
MHC_gene = loc[which(loc$chr == 6 & loc$start < 33448354 & loc$end > 28477797),]$'feature'
res = res[-which(res$gene %in% MHC_gene),]  # 10178

# specific gene set
specific_genes = res[order(res$'p.free.V'),'gene'][1:200]
write.table(specific_genes,'specific_set.txt',row.names = FALSE, col.names = FALSE, quote = FALSE)

# shared gene set
shared_genes = res[order(res$'p.free.hom_g2'),'gene'][1:200]
write.table(shared_genes,'shared_set.txt',row.names = FALSE, col.names = FALSE, quote = FALSE)

# GCTA gene set
gcta_res = fread('./op.greml')
gcta_res = gcta_res[-which(gcta_res$gene %in% MHC_gene),]
gcta_genes = gcta_res[order(gcta_res$p),'gene'][1:200]
write.table(gcta_genes,'GCTA_set.txt',row.names = FALSE, col.names = FALSE, quote = FALSE)

# DEG gene set
res_beta = res[,c('gene',paste0('free.ct_beta_',c('BIN','BMem','CD4ET','CD4NC','CD8ET','CD8NC','NK')))]
res_deg = as.data.frame(cbind(res_beta[,'gene'],apply(res_beta[,2:8],1,var)))
colnames(res_deg) = c('gene','deg')
res_deg$deg = as.numeric(res_deg$deg)
DEG_genes = res_deg[order(res_deg$deg,decreasing = TRUE),'gene'][1:200]
write.table(DEG_genes,'DEG_set.txt',row.names = FALSE, col.names = FALSE, quote = FALSE)



