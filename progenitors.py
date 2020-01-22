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


# ===================================================
#           Plotting
# ===================================================
def plot_ye(prog, ye, vline=None):
    fig, ax = plt.subplots(figsize=[8, 6])

    ax.plot(prog['radius'], prog['Ye'], label='Provided')
    ax.plot(prog['radius'], ye, label='Calculated')
    add_vline(ax, vline=vline, plot_type='ye')

    ax.set_xscale('log')
    ax.set_ylabel('Ye')
    ax.set_xlabel('radius (cm)')
    ax.legend()

    return fig, ax


def plot_abar(prog, abar, vline=None):
    fig, ax = plt.subplots(figsize=[8, 6])

    ax.plot(prog['radius'], prog['Abar'], label='Provided')
    ax.plot(prog['radius'], abar, label='Calculated')
    add_vline(ax, vline=vline, plot_type='abar')

    ax.set_xscale('log')
    ax.set_ylabel('Abar')
    ax.set_xlabel('radius (cm)')
    ax.legend()

    return fig, ax


def plot_x(prog, net, cr56=None, vline=None):
    fig, ax = plt.subplots(figsize=[8, 6])

    for row in net.itertuples():
        iso = row.isotope.capitalize()
        if iso == 'N':
            iso = 'neutrons'
        ax.plot(prog['radius']/1e5, prog[iso], label=iso)

    if cr56 is not None:
        ax.plot(prog['radius']/1e5, cr56, label='cr56')

    add_vline(ax, vline=vline, plot_type='x')

    ax.set_xscale('log')
    ax.set_ylabel('X')
    ax.set_xlabel('radius (km)')
    ax.legend()

    return fig, ax


def add_vline(ax, vline, plot_type):
    if vline is not None:
        ylims = None

        if plot_type in ('x', 'ye'):
            ylims = [0, 1]
        elif plot_type in ('abar',):
            ylims = [0, 60]

        ax.vlines(vline, ylims[0], ylims[1], ls='--', color='k')

