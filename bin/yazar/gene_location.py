import gzip, re
import numpy as np, pandas as pd

def read_gencode(f):
    res = {'gene':[], 'chr':[], 'start':[], 'end':[]}
    for line in gzip.open(f, 'rt'):
        if line[0] != '#':
            line = line.strip().split()
            if line[2] == 'gene' and re.search('chr',line[0]) and line[0][3:].isdigit():
                chr, start, end, info = int(line[0][3:]), int(line[3]), int(line[4]), line[-1]
                info = info.split(';')
                gene = info[0].split('.')[0][3:]
                res['gene'].append( gene ) 
                res['chr'].append( chr )
                res['start'].append(start) 
                res['end'].append(end)
    return( pd.DataFrame(res) )

def main():
    #
    genes = pd.read_table(snakemake.input.genes)
    gff = read_gencode( snakemake.input.gff )

    # merge
    genes = genes.merge(gff, left_on='feature', right_on='gene')
    genes = genes.drop('gene', axis=1)

    #
    genes.to_csv(snakemake.output.genes, sep='\t', index=False)

if __name__ == '__main__':
    main()

