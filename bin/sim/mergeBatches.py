import numpy as np, pandas as pd

def merge_dicts(reps, out):
    for key, value in reps[0].items():
        if isinstance(value, dict):
            out[key] = {}
            merge_dicts([rep[key] for rep in reps], out[key])
        else:
            out[key] = [rep[key] for rep in reps]

## list to array
def list2array(dic):
    for key, value in dic.items():
        #print(key)
        if key == 'out':
            dic[key] = np.array([])
            continue
        if isinstance(value, dict):
            list2array(value)
        elif isinstance(value, list):
            #print( value )
            dic[key] = np.array(value)

def main():
    # 
    he_outs = []
    reml_outs = []
    for f in snakemake.input.he:
        for line in open(f):
            out = np.load( line.strip(), allow_pickle=True ).item()
            he_outs.append( out )

    for f in snakemake.input.reml:
        for line in open(f):
            out = np.load( line.strip(), allow_pickle=True ).item()
            reml_outs.append( out )

    # merge HE and REML
    outs = []
    for he, reml in zip(he_outs, reml_outs):
        out = {'he': he, 'reml': reml}
        outs.append( out )

    # merge batches
    out = {}
    merge_dicts(outs, out)
    list2array( out )

    np.save( snakemake.output.out, out )

if __name__ == '__main__':
    main()
