import sys, os, re, tempfile, shutil
import numpy as np, pandas as pd
from cigma import util

def main():
    # read
    P = pd.read_table(sys.argv[1], index_col=0)
    op = pd.read_table(sys.argv[2], index_col=0)
    op.insert(loc=0, column='IID', value=op.index.to_numpy())
    grm = pd.read_table(sys.argv[3])
    op_pca = pd.read_table( sys.argv[4], index_col=0)
    op_pca = op_pca[[f'PC{i}' for i in range(1, int(sys.argv[5]) + 1)]]
    geno_pca = pd.read_table(sys.argv[6], index_col=0)
    geno_pca = geno_pca[[f'PC{i}' for i in range(1, int(sys.argv[7]) + 1)]]
    meta = pd.read_table(sys.argv[8], usecols=['individual','sex','age','pool'])
    meta = meta.drop_duplicates()
    meta.set_index('individual', inplace=True)
    snp_threshold = int(sys.argv[9])
    gene_out = sys.argv[10]
    out_f = sys.argv[11]

    # covariates
    tmpdir = tempfile.mkdtemp()
    tmp_f = os.path.join(tmpdir, 'tmp')
    covar_f = tmp_f + '.covar'
    qcovar_f = tmp_f + '.qcovar'
    pheno_f = tmp_f + '.pheno'
    batch_f = tmp_f + '.batch'
    mgrm_f = tmp_f + '.mgrm'

    # quantitative covariaes
    qcovar = op_pca.merge(geno_pca, left_index=True, right_index=True)
    qcovar = qcovar.merge(P.iloc[:, :-1], left_index=True, right_index=True)
    qcovar.insert(loc=0, column='IID', value=qcovar.index.to_numpy())
    qcovar.loc[op.index].to_csv(qcovar_f, header=False, sep='\t')

    # category covariates
    covar = meta[['sex']]
    # covar = meta[['sex', 'pool']] # correct batch as fixed effects
    covar = covar.merge(util.age_group(meta['age']), left_index=True, 
                        right_index=True)
    
    covar.insert(loc=0, column='IID', value=covar.index.to_numpy())
    covar.loc[op.index].to_csv(covar_f, header=False, sep='\t')

    # NOTE: make a fake GRM for batch effect
    # batch_grm = util.design(meta.index.to_numpy(), 
    #                         cat=meta['pool'].astype('category'), 
    #                         drop_first=False)
    # batch_grm = np.dot(batch_grm, batch_grm.T)  # NOTE: test if a coefficient would change V(G2)
    # batch_grm = pd.DataFrame(batch_grm, index=meta.index.to_numpy(), columns=meta.index.to_numpy())
    # batch_grm = batch_grm.where(np.tril(np.ones(batch_grm.shape)).astype(bool))
    # batch_grm = batch_grm.stack().reset_index()
    # batch_grm.columns = ['IID1', 'IID2', 'r']
    # batch_grm['nsnp'] = 1
    # find rank of each individual in meta
    # ranks = []
    # for ind in batch_grm['IID1']:
    #     x = np.argwhere(meta.index.to_numpy() == ind)
    #     if x.shape[0] != 1:
    #         sys.exit('wrong individual name')
    #     else:
    #         ranks.append(x[0][0] + 1)
    # batch_grm['IID1'] = ranks

    # ranks = []
    # for ind in batch_grm['IID2']:
    #     x = np.argwhere(meta.index.to_numpy() == ind)
    #     if x.shape[0] != 1:
    #         sys.exit('wrong individual name')
    #     else:
    #         ranks.append(x[0][0] + 1)
    # batch_grm['IID2'] = ranks
    # print(batch_grm.head())

    # batch_grm[['IID1', 'IID2', 'nsnp', 'r']].to_csv(batch_f + '.tmp.grm.gz', header=False, index=False, sep='\t')
    # tmp_meta = meta.reset_index()
    # tmp_meta[['individual', 'individual']].to_csv(batch_f + '.tmp.grm.id', index=False, header=False, sep='\t')
    # util.subprocess_popen(['gcta', '--grm-gz', batch_f + '.tmp', '--make-grm', '--out', batch_f])



    out = open(out_f, 'w')
    rng = np.random.default_rng(123)
    for index, row in grm.iterrows():
        gene, grm_f, snps = row['gene'], row['K'], row['snps']
        grm_prefix = os.path.splitext( os.path.splitext(grm_f)[0] )[0]

        if gene in op.columns:
            op[['IID',gene]].to_csv(pheno_f, header=False, sep='\t')
            gene_out_f = re.sub(r'/rep/',f'/{gene}/',gene_out)
        else:
            continue

        if snps > snp_threshold:
            os.makedirs(os.path.dirname(gene_out_f), exist_ok=True)
            # gcta
            # make multi grm file
            # with open(mgrm_f, 'w') as f:
                # f.write(f'{grm_prefix}\n{batch_f}')

            cmd1 = ['gcta', '--reml', '--reml-no-constrain', '--reml-maxit', '10000',
                    '--mgrm', mgrm_f, '--reml-no-lrt',
                    '--pheno', pheno_f, '--covar', covar_f, '--qcovar', qcovar_f, 
                    '--out', os.path.splitext(gene_out_f)[0]]
                    # '--reml-bendV', # doesn't not work
                    # '--grm', grm_prefix, 
            # reml-no-constrain doesn't work for reml-alg 2
            # reml-alg 1 can't converge with 10000 iters
            # reml-no-lrt reduced model has non-invertiable matrices
            cmd2 = ['gcta', '--reml', '--reml-no-constrain', 
                    '--reml-alg', '2',  '--reml-maxit', '10000',
                    '--grm', grm_prefix, 
                    '--pheno', pheno_f, '--covar', covar_f, '--qcovar', qcovar_f, 
                    '--out', os.path.splitext(gene_out_f)[0]]  
                    # '--mgrm', mgrm_f,
            # try:
            util.subprocess_popen( cmd2 )
            # except:
            #     sig = 0
            #     for i in range(10):
            #         try:
            #             sig = 0
            #             x = rng.random()
            #             y = 1 - x
            #             util.subprocess_popen( cmd1 + ['--reml-priors', str(x), str(y)])
            #             break
            #         except:
            #             sig = 1
            #     if sig == 1:
            #         continue

            out.write( gene_out_f + '\n' )
    
    # remove tmp files
    shutil.rmtree(tmpdir)


if __name__ == '__main__':
    main()
