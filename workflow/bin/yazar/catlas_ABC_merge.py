import numpy as np
import pandas as pd


def main():
    data = []
    
    for ct, f in zip(['B', 'T', 'NK', 'CD8', 'CD4'], 
                   [snakemake.input.B, snakemake.input.T, snakemake.input.NK, snakemake.input.CD8, snakemake.input.CD4]):
        df = pd.read_csv(f, sep='\t')
        df['ct'] = ct
        data.append(df)

    data = pd.concat(data, axis=0)

    # process
    data['chr'] = data['cCRE'].str.split(':').str[0].str.lstrip('chr')
    data['start'] = data['cCRE'].str.split(':').str[1].str.split('-').str[0].astype(int)
    data['end'] = data['cCRE'].str.split(':').str[1].str.split('-').str[1].astype(int)
    data = data.loc[data['chr'].isin([str(x) for x in range(1,23)])]

    data['promoter_chr'] = data['Promoter'].str.split(':').str[0].str.lstrip('chr')
    data['promoter_start'] = data['Promoter'].str.split(':').str[1].str.split('-').str[0].astype(int)
    data['promoter_end'] = data['Promoter'].str.split(':').str[1].str.split('-').str[1].astype(int)
    assert all(data['chr'] == data['promoter_chr'])


    data['gene'] = data['Gene Name'].str.split(':').str[0]
    data['transcript'] = data['Gene Name'].str.split(':').str[1]

    data = data[['cCRE', 'chr', 'start', 'end', 'promoter_start', 'promoter_end', 'gene', 'ct', 'ABC_score', 'Distance']]
    data = data.sort_values(['chr', 'start', 'end', 'ct']).reset_index(drop=True)

    # save
    data.to_csv(snakemake.output.abc, sep='\t', index=False)


if __name__ == '__main__':
    main()