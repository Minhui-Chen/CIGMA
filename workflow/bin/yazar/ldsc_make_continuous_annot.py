import numpy as np
import pandas as pd
from cigma import util
from pybedtools import BedTool


def main():
    # par
    out = np.load(snakemake.input.out, allow_pickle=True).item()
    location = pd.read_csv(snakemake.input.location, sep='\t')
    gcta = pd.read_table(snakemake.input.gcta, )
    bim = pd.read_table(snakemake.input.bim, names=['chr', 'rsid', 'cm', 'pos', 'a1', 'a2'])
    annot = bim[['chr', 'pos', 'rsid', 'cm']].copy().rename(columns={'chr':'CHR', 'pos':'BP', 'rsid':'SNP', 'cm':'CM'})
    window = int(float(snakemake.wildcards['window']))
    chr = int(snakemake.wildcards['chr'])

    # data process
    data = util.read_out(out, ['gene', 'free:hom_g2', 'free:v'])
    data = data.rename(columns={'free:hom_g2': 'hom_g2', 'free:v': 'v'})
    data['var_beta'] = out['free']['ct_beta'].var(axis=1)
    data['g'] = data['hom_g2'] + data['v']
    data = data.loc[data['g'] > 0] # remove negative genetic variance
    print(data.shape[0], 'genes with positive genetic variance')
    data['specificity'] = data['v'] / data['g']
    data = data.merge(gcta[['gene', 'h2']])

    data = data.merge(location[['feature', 'chr', 'start', 'end']], left_on='gene', right_on='feature').drop(columns=['feature'])
    data = data.loc[data['chr'] == chr]

    # remove MHC
    if chr == int(snkaemake.params.chr):
        mhc_start = int(snakemake.params.start)
        mhc_end = int(snakemake.params.end)
        data = data.loc[(data['start'] > mhc_end) | (data['end'] < mhc_start)]

    # add window
    data['raw_start'] = data['start']
    data['raw_end'] = data['end']
    data['start'] = np.maximum(1, data['start'] - window)
    data['end'] = data['end'] + window

    # make annotation
    vals = []
    for pos in annot["BP"]:
        nearby = data[
            (data["start"] <= pos) &
            (data["end"] >= pos)
        ]
        if not nearby.empty:
            # if AGG_METHOD == "max":
                # val = nearby["value"].max()
            # elif AGG_METHOD == "mean":
            val = nearby["specificity"].mean()
            # elif AGG_METHOD == "sum":
            #     val = nearby["value"].sum()
        else:
            val = 0.0
        vals.append(val)

    annot['specificity'] = vals
    annot.to_csv(snakemake.output.var, sep='\t', index=False, header=True)


if __name__ == '__main__':
    main()