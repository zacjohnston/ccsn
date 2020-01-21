import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt

prog_path = '/Users/zac/projects/data/progenitors/sukhbold_2018/mdotone'

def load_prog(mass, path=None, skiprows=3):
    """Load progenitor model

    mass : flt
    path : str
    skiprows : int
    """
    filepath = prog_filepath(mass, path=path)
    prog = pd.read_csv(filepath, skiprows=skiprows, delim_whitespace=True)
    return prog


def prog_filename(mass):
    """Return progenitor filename for given mass

    mass : flt
    """
    return f'{mass:.2f}.dat'


def prog_filepath(mass, path=None):
    """Return progenitor filepath for given mass

    mass : flt
    path : str
    """
    if path is None:
        path = prog_path

    filename = prog_filename(mass)
    filepath = os.path.join(path, filename)

    return filepath