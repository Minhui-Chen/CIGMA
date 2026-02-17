# src/cigma/datasets/__init__.py
"""
Dataset loading utilities for cigma package.
"""
import numpy as np
from importlib.resources import files
from pathlib import Path

def load_sim_data():
    """
    Load simulation dataset.
    
    Returns
    -------
    np.ndarray
        Simulation data array
        
    Examples
    --------
    >>> from cigma.datasets import load_sim_data
    >>> data = load_sim_data()
    >>> print(data.shape)
    """
    data_file = files('cigma.datasets').joinpath('sim.npy')
    with data_file.open('rb') as f:
        return np.load(f, allow_pickle=True).item()



# Convenience: expose at package level
__all__ = ['load_sim_data']