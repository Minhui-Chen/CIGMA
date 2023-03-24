snakemake@source('ctg_REML.source.R')
library(reticulate) # to use python
use_python('/apps/software/gcc-6.2.0/python/3.6.0/bin/python3')
np <- import('numpy')
os <- import('os')
source_python('bin/mystats.py')

# files
out_fs_name <- paste( snakemake@output[['out']], '.tmp', sep='' )

sink( out_fs_name )
for (i in snakemake@params[['batches']]) {
    cat( sub('/rep/', paste('/rep', as.character(i), '/', sep=''), snakemake@params[['out']]) )
    cat( '\n' )
}
sink()

cty_fs <- readLines( snakemake@input[['cty']] )
Z_fs <- readLines( snakemake@input[['Z']] )
ctvs_fs <- readLines( snakemake@input[['ctnu']] )
out_fs <- readLines( out_fs_name )

for ( i in 1:length(cty_fs) ) {
    cty <- as.matrix( read.table(cty_fs[i]) )
    Z <- as.matrix( read.table(Z_fs[i]) )
    ctvs <- as.matrix( read.table(ctvs_fs[i]) )

    hom <- screml_hom(cty, Z, ctvs)
    print(hom)
    print( system.time( free <- screml_free(cty, Z, ctvs) ) )
    print(free)
    #full <- screml_full(cty, Z, ctvs)
    C <- nrow( free[['V']] )
    lrt_p <- list( free_hom = lrt(free[['l']], hom[['l']], 2 * C)
                   #full_free = lrt( full[['l']], free[['l']], 2 * ( ((C-1)*C/2)-1 ) )
                  )

    # save
    py$out <- list(hom=hom, free=free)
    py$out[['lrt']] <- lrt_p
    dir.create( os$path$dirname(out_fs[i]), showWarnings = F, recursive = T )
    np$save( out_fs[i], py$out )
}

file.rename( out_fs_name, snakemake@output[['out']] )
