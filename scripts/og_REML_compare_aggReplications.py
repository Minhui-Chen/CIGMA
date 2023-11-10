import numpy as np
# waldNlrt out
def append_dicts(dic, dics):
    for key, value in dic.items():
        if isinstance(value, dict):
            if key not in dics.keys():
                dics[key] = {}
            append_dicts(value, dics[key])
        else:
            if key not in dics.keys():
                dics[key] = [value]
            else:
                dics[key].append(value)

out = {}
## aggregate
for f in snakemake.input.out:
    for line in open(f):
        npy = np.load(line.strip(), allow_pickle=True).item()
        append_dicts(npy, out)

## list to array
def list2array(dic):
    for key, value in dic.items():
        if isinstance(value, dict):
            list2array(value)
        elif isinstance(value, list):
            dic[key] = np.array(value)

list2array(out)
np.save(snakemake.output.out, out)
