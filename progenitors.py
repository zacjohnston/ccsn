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
    prog.rename(columns={'neutrons': 'Neutrons'}, inplace=True)  # consistent capitalise
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
def plot_ye(prog, ye=None, vline=None):
    fig, ax = plt.subplots(figsize=[8, 6])

    ax.plot(prog['radius'], prog['Ye'], label='Provided')
    add_calculated(ax=ax, vals=ye, prog=prog)
    add_vline(ax, vline=vline, plot_type='ye')

    ax.set_xscale('log')
    ax.set_ylabel('Ye')
    ax.set_xlabel('radius (cm)')
    ax.legend()

    return fig, ax


def plot_abar(prog, abar=None, vline=None):
    fig, ax = plt.subplots(figsize=[8, 6])

    ax.plot(prog['radius'], prog['Abar'], label='Provided')
    add_calculated(ax=ax, vals=abar, prog=prog)
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
        ax.plot(prog['radius'], prog[iso], label=iso)

    add_calculated(ax=ax, vals=cr56, prog=prog, label='cr56')
    add_vline(ax, vline=vline, plot_type='x')

    ax.set_xscale('log')
    ax.set_ylabel('X')
    ax.set_xlabel('radius (km)')
    ax.legend()

    return fig, ax


# ================================================================
#       Abundance Mapping
# ================================================================
# def map_abu(prog, net_0):
#     """
#     net_0 : table of isotopes *not* being mapped
#     """
#     abu = {}
#
#     abu['']
#
#     abu['fe56'] = prog['Fe56']
#
#     abu['cr56'] = 7 * (1 - sumx_0) - 14 * (prog['Ye'] - ye_0) - 0.5 * abu['fe56']
#
#     abu['ni56'] = 1 - sumx_0 - abu['fe56'] - abu['cr56']
#
#     return abu


def get_sums(prog, net):
    """net : table of isotopes to sum over
    """
    out = {}
    n_zones = len(prog)
    for key in ['sumx', 'ye', 'sumy']:
        out[key] = np.zeros(n_zones)

    for row in net.itertuples():
        iso = row.isotope.capitalize()
        x_i = np.array(prog[iso])

        out['sumx'] += x_i
        out['ye'] += x_i * (row.Z / row.A)
        out['sumy'] += x_i / row.A

    return out


# ================================================================
#       Convenience
# ================================================================
def add_calculated(ax, vals, prog, x_var='radius', label='Calculated'):
    if vals is not None:
        ax.plot(prog[x_var], vals, label=label)


def add_vline(ax, vline, plot_type):
    if vline is not None:
        ylims = None

        if plot_type in ('x', 'ye'):
            ylims = [0, 1]
        elif plot_type in ('abar',):
            ylims = [0, 60]

        ax.vlines(vline, ylims[0], ylims[1], ls='--', color='k')

