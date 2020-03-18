import matplotlib.pyplot as plt
import os
import sys

from nucleosynth.tracers import tracer


def plot(t, column, timesteps):
    path = '/home/zac/projects/codes/nucleosynth/temp/anim'
    fig, ax = plt.subplots(2, figsize=[10, 10])

    for i, timestep in enumerate(timesteps):
        print(f'\r{timestep}/{timesteps[-1]}', end='')
        ax[0].clear()
        ax[1].clear()
        t.plot_column(column, ax=ax[0], ylims=[4e8, 2e10], xlims=[0, 5])
        t.plot_sums(timestep, table='mass_frac', group='Z', ax=ax[1], title=False)

        time = t.columns['time'][timestep]
        ax[0].vlines(time, ymin=1e8, ymax=1e11, linewidth=1)

        filename = f'nuc_{i:04d}.png'
        filepath = os.path.join(path, filename)
        fig.savefig(filepath)

    print()
