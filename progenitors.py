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
def map_abu(prog, net_0, sums_0=None):
    """
    net_0 : table of isotopes *not* being mapped
    """
    abu = {}
    sums_0 = check_sums(sums_0, prog=prog, net=net_0)

    abu['fe56'] = prog['Fe56']

    abu['cr56'] = 7 * (1 - sums_0['sumx']) - 14 * (prog['Ye'] - sums_0['ye']) - 0.5 * abu['fe56']

    abu['ni56'] = 1 - sums_0['sumx'] - abu['fe56'] - abu['cr56']

    for key, val in abu.items():
        abu[key] = np.array(val)

    return abu


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


def plot_mapped(prog, net_0):
    fig, ax = plt.subplots(3, 1, figsize=[8, 6], sharex=True)

    for i in range(3):
        pass


def plot_mapped_sumx(prog, net_0, net, mapped_abu=None, ax=None, vline=None,
                     hline=None, x_var=None, sums_0=None, sums=None, xlims=None,
                     legend=True):
    x_var = check_xvar(x_var)
    sums_0 = check_sums(sums_0, prog=prog, net=net_0)
    sums = check_sums(sums, prog=prog, net=net)
    mapped_abu = check_mapped_abu(mapped_abu, prog=prog, net_0=net_0, sums_0=sums_0)

    ax = check_ax(ax)

    for iso in ['ni56', 'fe56', 'cr56']:
        ax.plot(prog[x_var], mapped_abu[iso], label=f'{iso} (mapped)')

    ax.plot(prog[x_var], sums['sumx'], ls='--', label='sumx (original)')
    ax.plot(prog[x_var], sums_0['sumx'], ls='--', label='sumx (partial)')

    sumx_final = sums_0['sumx'] + mapped_abu['ni56'] + mapped_abu['fe56'] + mapped_abu['cr56']
    ax.plot(prog[x_var], sumx_final, ls='--', label='sumx (final)')

    add_vline(ax, vline=vline, plot_type='x')
    add_hline(ax, hline=hline, prog=prog, x_var=x_var)
    config_ax(ax, legend=legend, xlims=xlims)
    return ax


def plot_mapped_ye(prog, net_0, net, mapped_abu=None, ax=None, vline=None,
                   hline=None, x_var=None, sums_0=None, sums=None, xlims=None,
                   legend=True):
    x_var = check_xvar(x_var)
    sums_0 = check_sums(sums_0, prog=prog, net=net_0)
    sums = check_sums(sums, prog=prog, net=net)
    mapped_abu = check_mapped_abu(mapped_abu, prog=prog, net_0=net_0, sums_0=sums_0)
    ax = check_ax(ax)

    ax.plot(prog[x_var], mapped_abu['ni56'] * 28/56, label=f'ni56 (mapped)')
    ax.plot(prog[x_var], mapped_abu['fe56'] * 26/56, label=f'fe56 (mapped)')
    ax.plot(prog[x_var], mapped_abu['cr56'] * 24/56, label=f'cr56 (mapped)')

    ax.plot(prog[x_var], prog['Ye'], ls='--', label='original')
    ax.plot(prog[x_var], sums_0['ye'], ls='--', label='partial')

    ye_final = sums_0['ye'] + mapped_abu['ni56'] * 28/56 + mapped_abu['fe56'] * 26/56 \
                 + mapped_abu['cr56'] * 24/56
    ax.plot(prog[x_var], ye_final, ls='--', label='final')

    add_vline(ax, vline=vline, plot_type='x')
    add_hline(ax, hline=hline, prog=prog, x_var=x_var)
    config_ax(ax, legend=legend, xlims=xlims, title='Ye')
    return ax


# ================================================================
#       Convenience
# ================================================================
def add_calculated(ax, vals, prog, x_var=None, label='Calculated'):
    x_var = check_xvar(x_var)
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


def add_hline(ax, hline, prog, x_var=None):
    x_var = check_xvar(x_var)

    if hline is not None:
        x = np.array(prog[x_var])
        xlims = [x[0], x[-1]]

        ax.hlines(hline, xlims[0], xlims[1], ls='--', color='k')


def check_xvar(x_var):
    if x_var is None:
        x_var = 'radius'
    return x_var


def check_sums(sums, prog, net):
    if sums is None:
        sums = get_sums(prog=prog, net=net)
    return sums


def check_mapped_abu(mapped_abu, prog, net_0, sums_0):
    if mapped_abu is None:
        mapped_abu = map_abu(prog=prog, net_0=net_0, sums_0=sums_0)
    return mapped_abu


def check_ax(ax):
    if ax is None:
        fig, ax = plt.subplots(figsize=[8, 6])
    return ax


def config_ax(ax, xscale='log', xlims=None, ylims=None, legend=True, title=None):
    ax.set_xscale(xscale)
    ax.set_title(title)
    xlims = check_xlims(xlims)
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)

    if legend:
        ax.legend()


def check_xlims(xlims):
    if xlims is None:
        xlims = [2e6, 1e10]
    return xlims
