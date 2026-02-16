import sys, shutil, os
import subprocess
import tempfile
import pandas as pd
import numpy as np


def extract_genotypes(bfile, chrom, start, end, maf, plink_path="plink"):
    """
    Extracts genotypes (0/1/2) from a VCF using PLINK for a given genomic window.

    Parameters
    ----------
    bfile : str
        Path to input binary PLINK file (without .bed/.bim/.fam extensions)
    chrom : str or int
        Chromosome number (e.g. '1')
    start : int
        Start position of window
    end : int
        End position of window
    maf : str
        Minor allele frequency threshold to filter SNPs
    plink_path : str, optional
        Path to plink executable (default: "plink")

    Returns
    -------
    geno_df : pandas.DataFrame
        DataFrame of genotype values (rows = samples, columns = SNPs)
    """
     # Create a temporary working directory
    tmpdir = tempfile.mkdtemp(prefix="plink_extract_")
    out_prefix = os.path.join(tmpdir, "subset")

    bim = pd.read_csv(bfile+'.bim', sep="\s+", header=None, names=["CHR","SNP","CM","BP","A1","A2"])
    subset = bim[(bim["CHR"] == chrom) & (bim["BP"].between(start, end))]
    if subset.shape[0] == 0:
        return np.empty((0,0))
    else:
        # Run PLINK to calculate MAF
        maf_out = os.path.join(tmpdir, "maf.frq")
        cmd_maf = [
            plink_path,
            "--bfile", bfile,
            "--chr", str(chrom),
            "--from-bp", str(start),
            "--to-bp", str(end),
            "--freq",
            "--out", os.path.join(tmpdir, "maf")
        ]
        print("Running:", " ".join(cmd_maf), flush=True)
        try:
            result = subprocess.run(cmd_maf, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print("STDOUT:", result.stdout.decode())
        except subprocess.CalledProcessError as e:
            print("Error occurred while running PLINK for MAF calculation")
            print("STDOUT:", e.stdout.decode())
            print("STDERR:", e.stderr.decode())
            sys.exit(1)
        maf_df = pd.read_csv(maf_out, sep="\s+")
        valid_snps = maf_df[maf_df["MAF"] > float(maf)]["SNP"].tolist()
        if len(valid_snps) == 0:
            shutil.rmtree(tmpdir)
            return np.empty((0,0))
        
        # Run PLINK to extract region and recode to additive (0/1/2)
        cmd = [
            plink_path,
            "--bfile", bfile,
            "--maf", str(maf),
            "--chr", str(chrom),
            "--from-bp", str(start),
            "--to-bp", str(end),
            "--recode", "A",  # Output additive coding (0/1/2)
            "--out", out_prefix
        ]
        print("Running:", " ".join(cmd), flush=True)
        try:
            result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print("STDOUT:", result.stdout.decode())
        except subprocess.CalledProcessError as e:
            print("Error occurred while running PLINK")
            print("STDOUT:", e.stdout.decode())
            print("STDERR:", e.stderr.decode())
            sys.exit(1)

        # Read the resulting .raw file
        raw_file = f"{out_prefix}.raw"
        df = pd.read_csv(raw_file, sep='\s+')

        # Drop PLINK metadata columns (FID, IID, etc.)
        geno_df = df.drop(columns=["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"])

        shutil.rmtree(tmpdir)

        return geno_df.values


def main():
    # par
    location_f = sys.argv[1]
    bed = sys.argv[2]
    bfile = bed.replace('.bed', '')
    out_f = sys.argv[3]
    chr = sys.argv[4]
    flank = int(sys.argv[5])
    maf = sys.argv[6]

    # read
    location_df = pd.read_table(location_f)
    location_df = location_df[location_df['chr'] == int(chr)]

    # extract genotypes
    data = {}
    for idx, row in location_df.iterrows():
        gene = row['feature']
        chrom = row['chr']
        start = max(1, int(row['start']) - flank)
        end = int(row['end']) + flank

        G = extract_genotypes(
            bfile,
            chrom,
            start,
            end,
            maf
        )

        data[gene] = {}
        data[gene]['G'] = G
        data[gene]['nsnps'] = G.shape[1]

    np.save(out_f, data)

if __name__ == '__main__':
    main()