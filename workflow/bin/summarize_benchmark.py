"""
Summarize runtime and memory of the HE benchmark runs (matrix-free vs dense he_ols).

Reads the merged benchmark outputs (one dict per file, as written by mergeBatches.py;
one entry per replicate in simulation or per gene in OneK1K):
    d['free']['time']   wall-clock seconds of free_HE
    d['peak_rss_mb']    peak resident memory of the job in MB (ru_maxrss)
and writes one CSV row per (analysis, implementation) with mean / SD / median of both.
"""
import numpy as np, pandas as pd


def fmt_time(s: float) -> str:
    """Seconds as in the paper's table, e.g. 0.02s, 4s, 20m, 2h."""
    if s < 1:
        return f'{s:.2g}s'
    if s < 60:
        return f'{s:.0f}s'
    if s < 3600:
        return f'{s / 60:.0f}m'
    return f'{s / 3600:.1f}h'


def fmt_mem(mb: float) -> str:
    """MB as GB, e.g. 0.2G, 17G."""
    gb = mb / 1024
    return f'{gb:.2g}G' if gb < 10 else f'{gb:.0f}G'


def summarize(row: dict, path: str) -> dict:
    d = np.load(path, allow_pickle=True).item()
    time = np.asarray(d['free']['time'], dtype=float)
    mem = np.asarray(d['peak_rss_mb'], dtype=float)
    out = dict(row)
    out.pop('file', None)
    out.update({
        'n': len(time),
        'time_mean_s': time.mean(), 'time_sd_s': time.std(ddof=1) if len(time) > 1 else np.nan,
        'time_median_s': np.median(time),
        'mem_mean_mb': mem.mean(), 'mem_sd_mb': mem.std(ddof=1) if len(mem) > 1 else np.nan,
        'mem_max_mb': mem.max(),
        'runtime': fmt_time(time.mean()), 'memory': fmt_mem(mem.mean()),
    })
    return out


def main():
    rows = snakemake.params.rows
    files = list(snakemake.input.files)
    assert len(rows) == len(files), 'rows and input files must correspond one-to-one'

    df = pd.DataFrame([summarize(row, f) for row, f in zip(rows, files)])
    df.to_csv(snakemake.output.csv, index=False)

    # human-readable copy in the log, in the layout of the paper's table
    cols = ['analysis', 'individuals', 'cell_types', 'implementation', 'n', 'runtime', 'memory']
    print(df[cols].to_string(index=False))


if __name__ == '__main__':
    main()
