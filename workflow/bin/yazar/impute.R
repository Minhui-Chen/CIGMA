source("bin/yazar/softImpute.R")

data <- read.table(snakemake@input[["data"]], sep='\t', header=F)
mat <- as.matrix(data)

out <- my_softImpute(mat, scale=T, seed=as.integer(snakemake@params[["seed"]]))
write.table(out$Y, snakemake@output[["data"]], sep='\t', quote=F, row.names=F, col.names=F)
