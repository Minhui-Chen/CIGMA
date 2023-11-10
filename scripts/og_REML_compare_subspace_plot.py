import numpy as np
import matplotlib.pyplot as plt


def main():
    #reml = np.load(snakemake.input.reml[0], allow_pickle=True).item()

    # plot
    fig, axes = plt.subplots(ncols=3, figsize=(15,4.8))

    #
    for reml_f in snakemake.input.reml:
        reml = np.load(reml_f, allow_pickle=True).item()

        for model in reml['reml'].keys():
            # hom2
            try:
                if model in ['null', 'hom', 'iid', 'free', 'full']:
                    axes[0].scatter(reml['reml'][model]['hom2'], reml['reml_proj'][model]['hom2'])
                else:
                    continue
            except:
                print(model,reml['reml'][model].keys())
                sys.exit()

            # V
            if 'V' in reml['reml'][model].keys():
                reml_V = np.array([np.diag(V) for V in reml['reml'][model]['V']]).flatten()
                reml_proj_V = np.array([np.diag(V) for V in reml['reml_proj'][model]['V']]).flatten()
                axes[1].scatter(reml_V, reml_proj_V)

            # W
            if 'W' in reml['reml'][model].keys():
                reml_W = np.array([np.diag(W) for W in reml['reml'][model]['W']]).flatten()
                reml_proj_W = np.array([np.diag(W) for W in reml['reml_proj'][model]['W']]).flatten()
                axes[2].scatter(reml_W, reml_proj_W)

    fig.savefig(snakemake.output.png)

if __name__ == '__main__':
    main()
