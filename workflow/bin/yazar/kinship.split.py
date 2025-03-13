import numpy as np
import gzip
import sys
import time

def main():
    batches = {}
    with open(snakemake.input.batch, 'r') as batch_file:
        for line in batch_file:
            chr, ctp_f = line.strip().split()
            batches[ctp_f] = int(chr)

    kinships = {}
    for chr, kinship_f in zip(snakemake.params.chrs, list(snakemake.input.kinship)):
        kinships[chr] = kinship_f

    k = 0
    for chr in snakemake.params.chrs:
        print(f"Processing chromosome {chr}", flush=True)
        data = np.load(kinships[chr], allow_pickle=True).item()
        chr_genes = np.array(data['gene'])

        for ctp_f, kinship_f in zip(snakemake.input.ctp, snakemake.output.kinship):
            # print(f"Time: {time.time()}", flush=True)
            # print(f"Processing {chr}", flush=True)
            if batches[ctp_f] != chr:
                continue
            
            with gzip.open(ctp_f, 'rt') as ctp_file:
                genes = np.array(ctp_file.readline().strip().split()[2:])
            
            kinship = {'ids': data['ids'], 'gene': [], 'K': [], 'nsnp': []}

            idx = np.isin(chr_genes, genes)
            if idx.sum() > len(genes) or idx.sum() > len(np.unique(chr_genes[idx])):
                print(genes)
                print(idx.sum())
                print(chr_genes[np.isin(chr_genes, genes)])
                print(genes[~np.isin(genes, chr_genes)])
                sys.exit('Missing genes')
            else:
                k += len(genes) - idx.sum()
                kinship['gene'] = chr_genes[idx]
                kinship['K'] = np.array(data['K'])[idx]
                kinship['nsnp'] = np.array(data['nsnp'])[idx]

            np.save(kinship_f, kinship)
            kinship = np.load(kinship_f, allow_pickle=True).item()
        
        # Explicitly delete large objects
        # del data, chr_genes, kinship
    
    print(f"Final k value: {k}", flush=True)

if __name__ == "__main__":
    main()
    print("Script completed successfully", flush=True)
    sys.stdout.flush()
    sys.stderr.flush()