import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

var_nu = pd.read_table(snakemake.input.var_nu, index_col=(0, 1))
nu = pd.read_table(snakemake.input.nu, index_col=(0, 1))

# intersection of ind-ct between var_nu and nu, restricting to pairs with # cell > threshold
shared_indices = var_nu.index.intersection(nu.index)
print(f'Before intersection: {var_nu.shape[0]}, {nu.shape[0]} pairs of ind-ct')
var_nu = var_nu.loc[shared_indices, nu.columns]
nu = nu.loc[shared_indices]

assert var_nu.index.equals(nu.index) 
assert var_nu.columns.equals(nu.columns)

print(f'After intersection: {var_nu.shape[0]} pairs of ind-ct')
print(f'After intersection: {var_nu.shape[1]} genes')

# calculate coefficient of variation
data = var_nu**(1/2) / nu

data.to_csv(snakemake.output.cv, sep='\t')

data = data.to_numpy().flatten()
# print(np.nanmax(data))
# print(np.nanmin(data))


fig, ax = plt.subplots()

# ax = axes[0]
# ax.hist(nu.to_numpy().flatten(), bins=100)
# ax.set_xlabel('nu')
# ax.set_ylabel('Frequency')


# ax = axes[1]
ax.hist(data, bins=30)
ax.set_xlabel('Coefficient of variation(nu)')
ax.set_ylabel('Frequency')

print(np.nanpercentile(data, (10, 50, 90)))
ax.axvline(np.nanpercentile(data, 10), ls='--', color='0.8')
ax.axvline(np.nanpercentile(data, 50), ls='--', color='0.8')
ax.axvline(np.nanpercentile(data, 90), ls='--', color='0.8')

# add coefficient of beta distribution
rng = np.random.default_rng()
mycolors = sns.color_palette()
ax.axvline(x=np.std(rng.choice([-1, 1], 10000) * rng.beta(.1, 1, 10000)),
                color=mycolors[1], ls='--', zorder=10, label='Beta(0.1, 1)')
ax.axvline(x=np.std(rng.choice([-1, 1], 10000) * rng.beta(.5, 1, 10000)),
                color=mycolors[2], ls='--', zorder=10, label='Beta(0.5, 1)')
# ax.axvline(x=np.std(rng.choice([-1, 1], 10000) * rng.beta(2, 1, 10000)),
                # color=mycolors[3], ls='--', zorder=10, label='Beta(2,5)')
# ax.axvline(x=np.std(rng.choice([-1, 1], 10000) * rng.beta(1, 1, 10000)),
                # color=mycolors[4], ls='--', zorder=10, label='Beta(2,3)')
ax.axvline(x=np.std(rng.choice([-1, 1], 1000) * rng.beta(20, 1, 1000)),
                color=mycolors[3], ls='--', zorder=9, label='Beta(20, 1)')
ax.legend()


fig.savefig(snakemake.output.png)
