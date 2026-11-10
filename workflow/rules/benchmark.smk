#########################################################################################
##   Summary of HE runtime and memory (matrix-free vs dense he_ols)
#########################################################################################
# Collects the merged benchmark outputs of
#   sim:     sim_HE_benchmark / sim_HE_dense_benchmark                 (Free model, base parameters)
#   OneK1K:  yazar_HE_free_[dense_]benchmark, yazar_trans_HE_free_[dense_]benchmark
# and writes one table of mean/SD runtime and memory per analysis and implementation.
#
#   snakemake results/benchmark/he_runtime_memory.csv

localrules: benchmark_he_summary


# simulation rows of the paper's table: (model, sample size, cell-type alpha)
benchmark_sim_rows = [
    ('free', '500', '2_2_2_2'),
    ('free', '1000', '2_2_2_2'),
    ('free', '1500', '2_2_2_2'),
    ('free', '2000', '2_2_2_2'),
    ('free21', '1000', '2_2'),
    ('free22', '1000', '2_2_2_2_2_2_2_2'),
]
benchmark_sim_vc = '0.13_0.13_0.11_0.04_0.08_0.51'  # base variance proportions of the Free model
benchmark_onek1k_inds, benchmark_onek1k_cts = 928, 7
benchmark_implementations = {'matrix-free': '', 'dense': '.dense'}


def _benchmark_rows():
    rows = []
    patterns = list(sim_paramspace.instance_patterns)  # same order as the rows of sim_params
    L, nL = config['sim']['L'], config['sim']['nL']
    for model, ss, a in benchmark_sim_rows:
        idx = sim_params.index[(sim_params['model'] == model) & (sim_params['ss'] == ss)
                               & (sim_params['a'] == a) & (sim_params['vc'] == benchmark_sim_vc)
                               & sim_params['V_diag'].str.startswith('1_1')]
        assert len(idx) == 1, f'expected one sim parameter row for {model}, {ss}, {a}'
        for impl, suffix in benchmark_implementations.items():
            rows.append({'analysis': f'Simulation ({model})', 'individuals': int(ss),
                         'cell_types': len(a.split('_')), 'implementation': impl,
                         'file': f'staging/sim/{model}/{patterns[idx[0]]}/L_{L}_nL_{nL}/out.he{suffix}.benchmark.npy'})

    # OneK1K, first (active) row of yazar.params.txt
    ypattern = next(iter(yazar_paramspace.instance_patterns))
    for analysis, folder in [('OneK1K cis', ''), ('OneK1K cis + trans', 'trans/')]:
        for impl, suffix in benchmark_implementations.items():
            rows.append({'analysis': analysis, 'individuals': benchmark_onek1k_inds,
                         'cell_types': benchmark_onek1k_cts, 'implementation': impl,
                         'file': f'staging/yazar/{ypattern}/{folder}he{suffix}.benchmark.npy'})
    return rows


benchmark_rows = _benchmark_rows()


rule benchmark_he_summary:
    input:
        files = [row['file'] for row in benchmark_rows],
    output:
        csv = 'results/benchmark/he_runtime_memory.csv',
    params:
        rows = benchmark_rows,
    script: '../bin/summarize_benchmark.py'
