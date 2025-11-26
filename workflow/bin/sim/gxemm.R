library(reticulate)
library(GxEMM)
use_condaenv("gxct")
np <- import("numpy")


# Load data
print("Loading data...")
data <- np$load(snakemake@input$data, allow_pickle=TRUE)  # data is a Python dict object
data <- data[[1]]  # extract the actual dict from the 0-d array

P <- as.data.frame(data[[names(data)[1]]][["P"]])
cts <- colnames(P)
C <- length(cts)

outs <- c(paste("gene", "sig2g_hom", paste0("sig2g_", 1:C, collapse='\t'), paste0("sig2e_", 1:C, collapse='\t'), "sig2e_hom", sep='\t'))

print(names(data))
for (name in names(data)) {
    y <- as.numeric(data[[name]][["y"]])
    P <- as.data.frame(data[[name]][["P"]])
    K <- as.matrix(data[[name]][["K"]])

    # Fit model
    out <- GxEMM_HE(y, P[,-1], K, P, Ztype='quant', gtype='free', etype='free') 
    # print(out)
    outs <- c(outs, paste(name, paste(out$sig2g, collapse='\t'), paste(out$sig2e, collapse='\t'), sep='\t'))

}


writeLines(outs, snakemake@output$out)
