.. These are examples of badges you might want to add to your README:
   please update the URLs accordingly

    .. image:: https://api.cirrus-ci.com/github/<USER>/CIGMA.svg?branch=main
        :alt: Built Status
        :target: https://cirrus-ci.com/github/<USER>/CIGMA
    .. image:: https://readthedocs.org/projects/CIGMA/badge/?version=latest
        :alt: ReadTheDocs
        :target: https://CIGMA.readthedocs.io/en/stable/
    .. image:: https://img.shields.io/coveralls/github/<USER>/CIGMA/main.svg
        :alt: Coveralls
        :target: https://coveralls.io/r/<USER>/CIGMA
    .. image:: https://img.shields.io/pypi/v/CIGMA.svg
        :alt: PyPI-Server
        :target: https://pypi.org/project/CIGMA/
    .. image:: https://img.shields.io/conda/vn/conda-forge/CIGMA.svg
        :alt: Conda-Forge
        :target: https://anaconda.org/conda-forge/CIGMA
    .. image:: https://pepy.tech/badge/CIGMA/month
        :alt: Monthly Downloads
        :target: https://pepy.tech/project/CIGMA
    .. image:: https://img.shields.io/twitter/url/http/shields.io.svg?style=social&label=Twitter
        :alt: Twitter
        :target: https://twitter.com/CIGMA

.. image:: https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold
    :alt: Project generated with PyScaffold
    :target: https://pyscaffold.org/

|

======
CIGMA
======

CIGMA is a Python package for decomposing cell type-shared and cell type-specific genetic effects on gene expression (eQTLs).
For a full description of the CIGMA model, please refer to the original paper: https://doi.org/??.??????/??????????.

This repository also contains the `Snakemake workflows <workflow/rules>`_ used for data analyses in our paper, including simulations and real-data applications.


Installation
============

Option 1: Conda environment (recommended)
------------------------------------------

The conda environment is defined in `environment.yml <envs/environment.yaml>`_, which includes Python 3.11.5 and all required dependencies:

.. code-block:: bash

   conda env create -n cigma -f envs/environment.yaml
   conda activate cigma

Option 2: pip install
---------------------

To install only the CIGMA Python package (should complete within seconds):

.. code-block:: bash

   pip install cigma

Python dependencies are listed in the ``install_requires`` section of `setup.cfg <setup.cfg>`_.

R dependencies (for CIGMA-REML)
-------------------------------

Running CIGMA-REML requires the following R packages:

.. code-block:: r

   install.packages(c("optparse", "numDeriv", "Matrix"))


Testing
=======

To verify the installation:

.. code-block:: bash

   python3 tests/test.py


Quick start
===========

CIGMA provides helper functions to generate model input from scRNA-seq data that have undergone quality control and normalization:

.. code-block:: python

   import scanpy as sc
   import pandas as pd
   from cigma import preprocess, fit

   # Read scRNA-seq data in h5ad format.
   # The data contains: a gene expression matrix (cells x genes) and cell metadata
   # with columns: 'cell' (cell IDs), 'ind' (individual IDs), 'ct' (cell types).
   ann = sc.read_h5ad(scRNA_h5ad_file)

   # Compute pseudobulk quantities:
   #   ctp  - cell type-pseudobulk matrix (individual-cell type pairs x genes)
   #   ctnu - cell-to-cell variation matrix (individual-cell type pairs x genes)
   #   P    - cell type proportion matrix (individuals x cell types)
   #   n    - cell count matrix (individuals x cell types)
   ctp, ctnu, P, n = preprocess.pseudobulk(
       ann=ann, ind_col='ind', ct_col='ct', cell_col='cell'
   )
   # or
   ctp, ctnu, P, n = preprocess.pseudobulk(
       X=ann.X, obs=ann.obs, var=ann.var, ind_col='ind', ct_col='ct', cell_col='cell'
   )

   # Remove genes or cell types to match the requirement of CIGMA that genes are expressed in all cell types.
   # For example, keep cell types CT1 and CT2 that all genes are expressed in, and remove other cell types.
   ctp = ctp.loc[ctp.index.get_level_values('ct').isin(['CT1', 'CT2'])]
   ctnu = ctnu.loc[ctnu.index.get_level_values('ct').isin(['CT1', 'CT2'])]
   P = P[['CT1', 'CT2']]
   # n = n.loc[:, ['CT1', 'CT2']]

   # Optional: scale so overall-pseudobulk has mean 0 and variance 1 across individuals.
   _, _, ctp, ctnu = preprocess.std(ctp, ctnu, P)

   # Fit a single gene with CIGMA-HE
   gene = 'GENE1'
   ctp_gene = ctp[gene].unstack().to_numpy()
   ctnu_gene = ctnu[gene].unstack().to_numpy()

   # Kinship matrix (individuals in same order as ctp/ctnu).
   # Can be computed from genotype data using GCTA or PLINK.
   K = pd.read_csv(K_file, index_col=0).to_numpy()

   out, pvalues = fit.free_HE(ctp_gene, K, ctnu_gene, P)


Input data
==========

See the `test script <tests/test.py>`_ for complete input data examples and usage.

.. _pyscaffold-notes:

Note
====

This project has been set up using PyScaffold 4.4. For details and usage
information on PyScaffold see https://pyscaffold.org/.
