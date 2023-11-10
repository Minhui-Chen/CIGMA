library(Matrix)
library(Seurat)
library(sctransform)
library(SeuratDisk)
library(reticulate)
use_condaenv('base')


# Set the warn option to display warnings
options(warn = 1)


# load python module
sc <- import("scanpy", convert=FALSE)
# sp <- import("scipy.sparse", convert=FALSE)


args <- commandArgs(trailingOnly = TRUE)



# AnnData to Seurat # not working, anndata to seurat has bug
# sceasy::convertFormat(args[1], from="anndata", to="seurat",
                    #    outFile=args[2])

# adata <- sc$read_h5ad(args[1], backed='r')
# adata <- sc$datasets$pbmc3k_processed()

# Instead, convert using seuratdisk
SeuratDisk::Convert(args[1], sub("h5ad", "h5seurat", args[1]), overwrite = TRUE)

# Create the Seurat object
# seurat <- CreateSeuratObject(counts = t(adata$X), meta.data = adata$obs)

seurat <- LoadH5Seurat(sub("h5ad", "h5seurat", args[1]))
# head(seurat@meta.data, 5)

# sctransform
seurat <- SCTransform(seurat, vst.flavor = "v2", vars.to.regress = c("percent.mt", "pool"), 
                        conserve.memory = TRUE, seed.use = 1448145)
saveRDS(seurat, file=args[2])

# Serat to AnnData
adata_seurat <- sc$AnnData(
    X   = t(GetAssayData(seurat)),
    obs = seurat[[]],
    var = GetAssay(seurat)[[]]
)

adata_seurat <- py_to_r(adata_seurat)

adata_seurat$write_h5ad(args[3])
