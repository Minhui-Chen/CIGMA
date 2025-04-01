.. These are examples of badges you might want to add to your README:
   please update the URLs accordingly

    .. image:: https://api.cirrus-ci.com/github/<USER>/GxCTMM.svg?branch=main
        :alt: Built Status
        :target: https://cirrus-ci.com/github/<USER>/GxCTMM
    .. image:: https://readthedocs.org/projects/GxCTMM/badge/?version=latest
        :alt: ReadTheDocs
        :target: https://GxCTMM.readthedocs.io/en/stable/
    .. image:: https://img.shields.io/coveralls/github/<USER>/GxCTMM/main.svg
        :alt: Coveralls
        :target: https://coveralls.io/r/<USER>/GxCTMM
    .. image:: https://img.shields.io/pypi/v/GxCTMM.svg
        :alt: PyPI-Server
        :target: https://pypi.org/project/GxCTMM/
    .. image:: https://img.shields.io/conda/vn/conda-forge/GxCTMM.svg
        :alt: Conda-Forge
        :target: https://anaconda.org/conda-forge/GxCTMM
    .. image:: https://pepy.tech/badge/GxCTMM/month
        :alt: Monthly Downloads
        :target: https://pepy.tech/project/GxCTMM
    .. image:: https://img.shields.io/twitter/url/http/shields.io.svg?style=social&label=Twitter
        :alt: Twitter
        :target: https://twitter.com/GxCTMM

.. image:: https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold
    :alt: Project generated with PyScaffold
    :target: https://pyscaffold.org/

|

======
GxCTMM
======


|GxCTMM| is a Python package for the decomposition of cell type-shared and -specific eQTLs using the CIGMA model.
For a full description of CIGMA, please refer to the original paper: https://doi.org/10.1101/2023.08.01.551679.

This repository contains scripts for data analyses in our paper. [Snakefiles](workflows/rules) contains steps for running CIGMA model on simulated and real data.

.. * Download GWAS data from ... and update the path in the [config](config/config.yaml) file.
.. * Download LDSC: git clone https://github.com/bulik/ldsc.git


Installation
======
The conda env is defined in the [environment.yml](environment.yml) file.
To create the conda environment, run:
```bash
conda env create -n gxct -f environment.yml
conda activate gxct
```
To only install the GxCTMM Python package, run:
```bash
pip install gxctmm
```
To run the tests, run:
```bash
python3 tests/test.py
```

.. _pyscaffold-notes:

Input data
======
Please check the [test script](tests/test.py) for CIGMA input data and running examples.

Note
====

This project has been set up using PyScaffold 4.4. For details and usage
information on PyScaffold see https://pyscaffold.org/.
