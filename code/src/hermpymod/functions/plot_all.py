import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import warnings 

from hermpymod.classes.panels import PlanarplotPanel


def plot_all_ephemeris(times, units="Mercury Radii", crossings=None ,encounters=None, plot=False, ax=None, color="C0", scatter=True, mercury=True, label=None, add_legend=True, frame="MSO"):

    if ax == None:
        ax_true = True 
        fig = plt.figure()
        gs = GridSpec(2, 3, figure=fig)

        ax = [
        fig.add_subplot(gs[0, 0]),
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[0, 2]), 
        fig.add_subplot(gs[1, :]),
        ]

    else:
        ax_true = False
    
    plane_plot = PlanarplotPanel(times, plane='All', units=units, crossings=crossings, encounters=encounters, color=color, scatter=scatter, mercury=mercury, label=label, frame=frame)

    plane_plot._plot_on(ax)

    handles, _ = ax[0].get_legend_handles_labels()

    for axis in ax:
        axis.xaxis.label.set_size(14)
        axis.yaxis.label.set_size(14)
        axis.tick_params(axis='both', labelsize=12)
        axis.legend().remove()

    if ax_true:
        fig.suptitle(f"MESSENGER SPICE data, taken from {times[0]}-{times[-1]}", fontsize=18)
        if add_legend:
            fig.legend(handles=handles,  bbox_to_anchor=(1.0,0.5))

    if not ax_true:
        return None
    elif plot:
        return fig, ax, plane_plot
    else:
        return fig, ax

