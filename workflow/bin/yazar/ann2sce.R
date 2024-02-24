library(reticulate)
library(SingleCellExperiment)
library(Matrix)
use_condaenv('base')

# load python module
sc <- import("scanpy", convert=FALSE)

args <- commandArgs(trailingOnly = TRUE)
input_f <- args[1]
output_f <- args[2] 

adata <- sc$read_h5ad(input_f, backed='r')
print(adata$obs$head())
print(adata$var$head())

# make SCE  
# SparseDataset workaround: use H5ADMatrix to read (auto tranposed matrix relative to python); otherwise, zellkonverter
sce <- SingleCellExperiment(list(counts=HDF5Array::H5ADMatrix(input_f)),  
                            colData=py_to_r(adata$obs),
                            rowData=py_to_r(adata$var))  

print(head(colData(sce)))
print(head(rowData(sce)))

saveRDS(sce, output_f)
