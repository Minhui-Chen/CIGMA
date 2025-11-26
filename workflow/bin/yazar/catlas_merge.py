import os
import numpy as np
import pandas as pd


def main():
    dirname = os.path.dirname(snakemake.input.bed)
    bed_fs = [os.path.join(dirname, f) for f in os.listdir(dirname) if f.endswith('.bed')]
    dfs = []
    for f in bed_fs:
        df = pd.read_table(f, names=['chr', 'start', 'end', 'cre', 'x', 'y'])
        df.drop(columns=['x', 'y'], inplace=True)
        df['cell_type'] = os.path.basename(f).replace('.bed', '')
        dfs.append(df)
    df = pd.concat(dfs, axis=0)
    df['chr'] = df['chr'].str.replace('chr', '')
    # remove sex chromosomes
    df = df[df['chr'].isin(np.arange(1, 23).astype(str))]
    df = df.sort_values(['chr', 'start', 'end', 'cre', 'cell_type'])
    df.to_csv(snakemake.output.cre, sep='\t', index=False, header=True)


if __name__ == '__main__':
    main()