import matplotlib.pyplot as plt
import os
import numpy as np
import sys

from nucleosynth.tracers import tracer


def plot(t, column, times):
    path = '/home/zac/projects/codes/nucleosynth/temp/anim'
    fig, ax = plt.subplots(2, figsize=[10, 10])
    n = len(times)

    for i, time in enumerate(times):
        print(f'\r{i}/{n}', end='')
        ax[0].clear()
        ax[1].clear()

        # idx = find_nearest_idx(t.columns['time'], time)
        idx = np.searchsorted(t.columns['time'], time)

        t.plot_column(column, ax=ax[0], ylims=[1.5e9, 1.2e10], xlims=[0.1, 1.5])
        t.plot_sums(idx, table='mass_frac', group='Z', ax=ax[1], title=False,
                    xlims=[-1, 40], ylims=[1e-35, 1e1])

        ax[0].vlines(time, ymin=1e8, ymax=1e11, linewidth=1)

        filename = f'nuc_{i:04d}.png'
        filepath = os.path.join(path, filename)
        fig.savefig(filepath)

    print()


def find_nearest_idx(array, value):
    """Return idx for the array element nearest to the given value

    Note: array assumed to be monotonically increasing (not enforced),
          will use the first element that exceeds the given value

    parameters
    ----------
    array : 1D array
        array to search
    value : float
        value to look for in array
    """
    idx = np.searchsorted(array, value)
    if np.abs(value - array[idx - 1]) < np.abs(value - array[idx]):
        return idx - 1
    else:
        return idx




