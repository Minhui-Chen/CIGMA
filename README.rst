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


CIMGA is a Python package for the decomposition of cell type-shared and -specific eQTLs using the CIGMA model.
For a full description of CIGMA, please refer to the original paper: https://doi.org/10.1101/2023.08.01.551679.

This repository contains scripts for data analyses in our paper. `Snakemake files <workflow/rules>`_ contain steps for running CIGMA model on simulated and real data.

.. * Download GWAS data from ... and update the path in the [config](config/config.yaml) file.
.. * Download LDSC: git clone https://github.com/bulik/ldsc.git


Installation
============
The conda env is defined in the `environment.yml <env/environment.yml>`_ file.

To create the conda environment, run:

```bash
conda env create -n cigma -f env/environment.yml
conda activate cigma
```


To only install the CIGMA Python package, run:

```bash
pip install cigma
```


To run the tests, run:

```bash
python3 tests/test.py
```


.. _pyscaffold-notes:

Input data
==========
Please check the `test script <tests/test.py>`_ for CIGMA input data and running examples.

Note
====

This project has been set up using PyScaffold 4.4. For details and usage
information on PyScaffold see https://pyscaffold.org/.
