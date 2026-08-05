import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.colors as mcolors

from hermpy.net import ClientSPICE, ClientMESSENGER
from hermpy.data import parse_messenger_mag
from hermpy.utils import Constants

from astropy.time import Time
from sunpy.time import TimeRange

import warnings 

from hermpymod.classes.panels import PlanarplotPanel


c = ClientMESSENGER()

Zd = Constants.DIPOLE_OFFSET


def plot_all_ephemeris(times, units="Mercury Radii", crossings=None ,encounters=None, plot=False, ax=None, color="C0", scatter=True, mercury=True, label=None, add_legend=True, frame="MSO", downsampled=True, resolution=None, boundaries=False, mag=False):

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
    
    plane_plot = PlanarplotPanel(times, plane='All', units=units, crossings=crossings, encounters=encounters, color=color, scatter=scatter, mercury=mercury, label=label, frame=frame, downsampled=downsampled, resolution=resolution, boundaries=boundaries)

    plane_plot._plot_on(ax)

    handles, _ = ax[0].get_legend_handles_labels()

    for axis in ax:
        axis.xaxis.label.set_size(14)
        axis.yaxis.label.set_size(14)
        axis.tick_params(axis='both', labelsize=12)
        axis.legend().remove()

    if mag:
        mag_time_range = TimeRange(Time(times[0]).to_datetime(), Time(times[-1]).to_datetime())
        c.query(mag_time_range, "MAG")
        mag_data_encounter = c.fetch()

        mag_table : QTable = parse_messenger_mag(mag_data_encounter, mag_time_range)

        slicing = len(mag_table)//100

        mag_table=mag_table[::slicing]

        mag_total = np.sqrt(mag_table["Bx"]**2 + mag_table["By"]**2 + mag_table["Bz"]**2).value
        
        pos_col = [col for col in mag_table.keys() if "X" in col or "Y" in col or "Z" in col]
        
        for col in pos_col:
            mag_table[col] = mag_table[col].to(units)

        if frame == "MSM":
            mag_table["Z MSO"] = mag_table["Z MSO"] - Zd.to(units)

        mag_table["Rho MSO"] = np.sqrt(mag_table["Y MSO"]**2 + mag_table["Z MSO"]**2)
        mag_table["Brho"] = np.sqrt(mag_table["By"]**2 + mag_table["Bz"]**2)

        mag_config = [
                ("X MSO", "Y MSO", "Bx", "By", ax[0]),
                ("X MSO", "Z MSO", "Bx", "By", ax[1]),
                ("Y MSO", "Z MSO", "Bx", "By", ax[2]),
                ("X MSO", "Rho MSO", "Bx", "Brho", ax[3]),
                ]

        norm = mcolors.Normalize(vmin=0, vmax=120)


        for xlab, ylab, bxlab, bylab, axis in mag_config:
            mag_cm = axis.quiver(mag_table[xlab].value, mag_table[ylab].value, mag_table[bxlab].value, mag_table[bylab].value, mag_total, cmap='viridis', scale_units='width', scale=850, norm=norm)


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

