import shutil
import sys, os, re, tempfile
import numpy as np, pandas as pd
from cigma import util


def check_convergence(log_f):
    lines = open(log_f, 'r').readlines()
    for i in range(len(lines)):
        if "Iter.\tlogL\tV(G)\tV(e)" in lines[i]:
            break  # Stop checking after iterations start
    if i == len(lines) - 1:
        sys.exit("Error: Iteration log not found in GCTA output.")
    k = 0
    for j in range(i+1, len(lines)):
        line = lines[j].strip().split()
        if line[0].isdigit():
            new_k = int(line[0])
            if new_k != k + 1:
                sys.exit(f"Error: Iteration numbers are not consecutive at line {j}: {line}")
            k = new_k
        elif "Log-likelihood ratio converged." in lines[j]:
            break
        else:
            sys.exit(f"Error: Unexpected line format at line {j}: {line}")
    return k


def main():
    #
    data = np.load(sys.argv[1], allow_pickle=True).item()
    outs = []

    tmpdir = tempfile.mkdtemp()
    tmp = os.path.join(tmpdir, 'tmp')
    for key in data.keys():
        C = data[key]['P'].shape[1]

        prefix = tmp + f'.{key}'
        pheno_f = prefix + '.pheno'

        # gcta command
        cmd1= ['gcta', '--reml', '--reml-no-constrain',
                '--grm-gz', prefix, '--pheno', pheno_f,
                '--out', prefix]
        # NOTE: reml alg
        cmd2 = ['gcta', '--reml', '--reml-no-constrain', 
                '--reml-alg', '2', '--reml-maxit', '10000',
                '--grm-gz', prefix, '--pheno', pheno_f,
                '--out', prefix]


        # generate gxta grm
        inds = util.grm_matrix2gcta_grm_gz(data[key]['K'], prefix, nsnp=data[key]['G'].shape[1])
        
        gremls = []
        for i in range(C):
            y = data[key]['Y'][:,i]
            pheno = np.column_stack((np.array([inds, inds]).T, y))
            np.savetxt(pheno_f, pheno, fmt='%s', delimiter='\t')

            # run GCTA
            util.subprocess_popen( cmd2 )
            k = check_convergence(prefix + '.log')

            greml = util.read_greml(prefix + '.hsq')
            greml['converged'] = "Log-likelihood ratio converged at iteration " + str(k)
            gremls.append(greml)
        
        # merge into one dict
        greml = {}
        util.merge_dicts(gremls, greml)

        # OP h2
        y = data[key]['y']
        pheno = np.column_stack((np.array([inds, inds]).T, y))
        np.savetxt(pheno_f, pheno, fmt='%s', delimiter='\t')
        qcovar_f = prefix + '.qcovar'
        P = data[key]['P']
        qcovar = np.column_stack((np.array([inds, inds]).T, P[:, :-1]))  # NOTE: gcta has a mean column
        np.savetxt(qcovar_f, qcovar, fmt='%s', delimiter='\t')
        util.subprocess_popen( cmd2 + ['--qcovar', qcovar_f] )
        greml_op = util.read_greml(prefix + '.hsq')
        k = check_convergence(prefix + '.log')
        greml_op['converged'] = "Log-likelihood ratio converged at iteration " + str(k)


        # save
        outs.append({'gene': key, 'greml_ctp': greml, 'greml_op': greml_op})

    shutil.rmtree(tmpdir)
    
    np.save(sys.argv[2], outs)


if __name__ == '__main__':
    main()
