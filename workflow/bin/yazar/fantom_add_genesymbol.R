library(clusterProfiler)
library(AnnotationHub)

# OrgDb
hub <- AnnotationHub()
q <- query(hub, "org.Hs.eg.db")
id <- q$ah_id[length(q)]
org.Hs.eg.db <- hub[[id]]

# read
df <- read.table(snakemake@input[['promoter']], header = TRUE)

# convert
gene <- bitr(na.omit(df[, 'entrezgene_id']), fromType="ENTREZID", toType=c("SYMBOL", "ENSEMBL"), OrgDb=org.Hs.eg.db, drop=F)
write.table(gene, snakemake@output[['genemeta']], sep='\t', quote=FALSE, row.names=FALSE)