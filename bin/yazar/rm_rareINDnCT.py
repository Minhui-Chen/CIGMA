import numpy as np, pandas as pd
from memory_profiler import profile

@profile
def main():
    # read
    ctp = pd.read_table(snakemake.input.ctp, dtype='float32', index_col=(0,1))
    print( ctp.head(), flush=True )
    ctnu = pd.read_table(snakemake.input.ctnu, dtype='float32', index_col=(0,1))
    print( ctnu.head(), flush=True )
    n = pd.read_table(snakemake.input.n, index_col=0)
    print( n.head(), flush=True )

    # 
    n_inds = n.shape[0]
    cts = n.columns.to_numpy()
    
    # find common cts with less than 50% missing inds
    print( '1', flush=True )
    common = n.median()[n.median() > int(snakemake.wildcards.ct_min_cellnum)].index.tolist()

    # filter common cts
    ctp = ctp.loc[ctp.index.get_level_values('ct').isin(common)] 
    ctnu = ctnu.loc[ctnu.index.get_level_values('ct').isin(common)] 
    n = n[common]

    # filter rare inds
    print( '1' )
    n_cells = n.sum(axis=1)
    common_ids = n_cells[n_cells > int(snakemake.wildcards.ind_min_cellnum)].index.tolist()
    ctp = ctp.loc[ctp.index.get_level_values('ind').isin(common_ids)]
    ctnu = ctnu.loc[ctnu.index.get_level_values('ind').isin(common_ids)]
    n = n[n.index.isin(common_ids)]

    # save n and P
    print( '1' )
    n.to_csv(snakemake.output.n, sep='\t')
    P = n.div(n.sum(axis=1), axis=0)
    P.to_csv(snakemake.output.P, sep='\t')

    # filter rare ind-cts
    print( '1' )
    n = n.stack()
    common_cts = n[n > int(snakemake.wildcards.ct_min_cellnum)].index
    ctp = ctp.loc[ctp.index.isin( common_cts )]
    ctnu = ctnu.loc[ctnu.index.isin( common_cts )]

    # 
    print( '1' )
    ctp.to_csv(snakemake.output.ctp, sep='\t')
    ctnu.to_csv(snakemake.output.ctnu, sep='\t')

if __name__ == '__main__':
    main()


