library(Matrix)
library(zellkonverter)
# library(reticulate)
library(SingleCellExperiment)
# use_condaenv('base')

# load python module
# sc <- import("scanpy", convert=FALSE)

args <- commandArgs(trailingOnly = TRUE)
input_rds <- args[1]
transform <- args[2] 
transformation <- strsplit(transform, '.', fixed=T)[[1]][1]
alpha <- paste(strsplit(transform, '.', fixed=T)[[1]][-1], collapse='.')
seed <- args[3]
output_f <- args[4]

set.seed(seed)

# get a list called `all_transformations` that contains 
# all transformations as ready to call functions
# Furthermore, it initializes the `make_knn_graph` function
source("transformGamPoi-Paper/benchmark/src/transformations/transformation_helper.R")  # TODO


######### Start Transformation #######
sce <- readRDS(input_rds)
print(head(colData(sce)))
print(head(rowData(sce)))
UMI <- as.matrix(counts(sce))

# remove cells / genes not expressed at all
print('1')
expressed_cells <- matrixStats::colSums2(UMI) > 0
expressed_genes <- matrixStats::rowSums2(UMI) > 0
sce <- sce[expressed_genes, expressed_cells]
UMI <- as.matrix(counts(sce))
print(sce)

if(alpha == "global"){
  alpha <- "global"
}else if(! is.na(suppressWarnings(readr::parse_double(alpha, na = character(0L))))){
  alpha <- readr::parse_double(alpha)
}else if(! is.na(suppressWarnings(readr::parse_logical(alpha, na = character(0L))))){
  alpha <- readr::parse_logical(alpha)
}else{
  stop("Cannot parse alpha=", alpha)
}

sf <- MatrixGenerics::colSums2(UMI)
sf <- sf / mean(sf)

print('4')
trans_UMI <- all_transformations[[transformation]](UMI, sf, alpha)
trans_UMI <- Matrix(trans_UMI)  # make sparse # TODO: sparse matrix format csc vs csr?

# trans to ann
print('5')
counts(sce) <- trans_UMI
saveRDS(trans_UMI, file=paste0(tools::file_path_sans_ext(output_f), ".rds"))
adata <- SCE2AnnData(sce)

print('6')

adata$write_h5ad(output_f)