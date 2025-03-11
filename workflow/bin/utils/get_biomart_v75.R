library(biomaRt)

d_archives <- listEnsemblArchives()

url <- d_archives[d_archives$version == '75', 'url']

mart <- useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', host=url)

chrs <- 1:22
bm <- getBM(c("ensembl_gene_id", "hgnc_symbol", "hgnc_id", "entrezgene", "chromosome_name", "start_position", "end_position", "transcript_start", "transcript_end", "strand"), 
            filters="chromosome_name", values=chrs, mart=mart)

# bm <- subset(bm, start_position == transcript_start | end_position == transcript_start)
# bm <- unique(bm)
# print(bm)
write.table(bm, snakemake@output[['gene_info']], sep='\t', quote=FALSE, row.names=FALSE)