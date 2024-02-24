import sys, os, time
import scanpy as sc

from ctmm import preprocess

import tracemalloc, linecache

def display_top(snapshot, key_type='lineno', limit=3):
    snapshot = snapshot.filter_traces((
        tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
        tracemalloc.Filter(False, "<unknown>"),
    ))
    top_stats = snapshot.statistics(key_type)

    print("Top %s lines" % limit)
    for index, stat in enumerate(top_stats[:limit], 1):
        frame = stat.traceback[0]
        # replace "/path/to/module/file.py" with "module/file.py"
        filename = os.sep.join(frame.filename.split(os.sep)[-2:])
        print("#%s: %s:%s: %.1f KiB"
              % (index, filename, frame.lineno, stat.size / 1024))
        line = linecache.getline(frame.filename, frame.lineno).strip()
        if line:
            print('    %s' % line)

    other = top_stats[limit:]
    if other:
        size = sum(stat.size for stat in other)
        print("%s other: %.1f KiB" % (len(other), size / 1024))
    total = sum(stat.size for stat in top_stats)
    print("Total allocated size: %.1f KiB" % (total / 1024))

h5ad = sys.argv[1]
genes = [line.strip() for line in open(sys.argv[2])]
path = sys.argv[3]

#tracemalloc.start()

data = sc.read_h5ad(h5ad, backed='r')
print('Finish reading', flush=True)

subdata = data[(~data.obs.ind_cov.isna()) & (~data.obs.author_cell_type.isna()), genes]
print('Subsetting', flush=True)

#if os.path.exists('test.h5ad'):
#    os.remove('test.h5ad')
#subdata.write('test.h5ad', compression='gzip')
#snapshot = tracemalloc.take_snapshot()
#display_top(snapshot)

ctp, ctnu, P = preprocess.pseudobulk(ann=subdata, ind_cut=100, ct_cut=10, 
        ind_col='ind_cov', ct_col='author_cell_type')
#ctp, ctnu, P = preprocess.pseudobulk(meta=subdata.obs.reset_index(drop=False,names='cell'), counts=subdata.to_df().T, 
#        matched=True, ind_cut=100, ct_cut=10, 
#        ind_col='ind_cov', ct_col='author_cell_type')
print( ctp.head() )
print('Done', flush=True)


ctp.to_csv(path+'.ctp.gz', sep='\t')
ctnu.to_csv(path+'.ctnu.gz', sep='\t')
P.to_csv(path+'.P.gz', sep='\t')
