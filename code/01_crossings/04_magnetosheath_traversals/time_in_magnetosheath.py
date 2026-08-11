import os
from pathlib import Path
os.chdir(Path(__file__).resolve().parent)

import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import numpy as np
from astropy.table import QTable, vstack, Column
from hermpy.plotting import Panel, TimeseriesPanel
from hermpy.net import ClientSPICE, ClientMESSENGER
from astropy.time import Time
import spiceypy as spice
from hermpy.utils import Constants as c
import datetime as dt
from matplotlib.lines import Line2D
from astropy.table import QTable, hstack, vstack
import sys
import os

from hermpymod.classes.panels import PlanarplotPanel,  HistogramPanel
from hermpymod.functions.plot_all import plot_all_ephemeris
from hermpymod.functions.ephemeris_downsampler import parse_crossing_list


home_dir = os.getenv('HOME')
from hermpymod.paths import DATA_DIR
data_dir = DATA_DIR
img_dir = "../../../plots_and_images/"

plt.style.use(img_dir + "presentation.mplstyle")

spice_client = ClientSPICE()

hollman_crossing_list = parse_crossing_list()
hollman_crossing_list["UTC"] = Time(hollman_crossing_list["UTC"]).to_datetime()


"""
Function to pair the crossings just before and after the magnetosheath
"""


def delta_t_magenetosheath(crossing_list):
    crossing_time = Time(crossing_list["UTC"]).to_datetime()
    crossing_label = crossing_list["Label"]

    crossings_delta_t = [(crossing_time[i+1] - crossing_time[i]) for i in range(len(crossing_time) -1)]

    crossings_delta_t = [i.total_seconds()/3600 for i in crossings_delta_t]

    ms_out = []
    ms_in = []


    for idx, crossing in enumerate(crossing_list):
        if crossing["Label"] == "MP_OUT":
            if crossing_label[idx + 1] == "BS_OUT":
                ms_out.append([crossing, crossing_list[idx + 1]])
            else:
                continue
        elif crossing["Label"] == "BS_IN":
            if crossing_label[idx + 1] == "MP_IN":
                ms_in.append([crossing, crossing_list[idx + 1]])
            else:
                continue
        

    return crossings_delta_t, ms_out , ms_in

configs = [
        (hollman_crossing_list, "Hollman")
        ]

# ~5 Minute bins
bins = np.arange(0,8,0.084)
for data, name in configs:
    dt, ms_out, ms_in = delta_t_magenetosheath(data)

    """
    Plots to make, 
    inbound magnetosheath traversals
    outbound magnetosheath traversals
    """

    data_sets = [
            (ms_in, r"BS IN \& MP IN", "purple", "inbound"),
            (ms_out, r"MP OUT \& BS OUT", "pink", "outbound"),
            ]
    
    plots = []

    id_start =0
    id_end = 100

    for data_type, label, color, direction in data_sets:

        # Magnetosheath traversal duration
        ms_dt = np.array([(i[-1]["UTC"] - i[0]["UTC"]).total_seconds()/3600 for i in data_type])

        # Take shortest 100 traversals
        shortest_indices = np.argsort(ms_dt)[id_start:id_end]

        ms_short = []


        for idx in shortest_indices:
            ms_short.append(data_type[idx])


        fig = plt.figure()
        gs = GridSpec(2, 3, figure=fig)

        ax = [
        fig.add_subplot(gs[0, 0]),
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[0, 2]), 
        fig.add_subplot(gs[1, :]),
        ]

        # Stack crossings for plotting
        ms_crossings = [vstack(i) for i in ms_short]
        ms_crossings = vstack(ms_crossings)
        

        # Plot shortest magnetosheath traversals
        with spice_client.KernelPool():
            for idx in shortest_indices:
                ms = data_type[idx]
                time_range = [ms[0]["UTC"].isoformat(),ms[-1]["UTC"].isoformat()]
                plot_all_ephemeris(time_range, color='k', ax=ax, scatter=False, downsampled=False, resolution=10, frame="MSM")

        # Plot crossings
        plot_all_ephemeris([data_type[0][0]["UTC"].isoformat(),data_type[-1][-1]["UTC"].isoformat()], color='none', ax=ax, crossings=ms_crossings, frame="MSM")

        fig.suptitle(label)

        for axis in ax[:-2]:
            axis.set_xlim(0,2)
            axis.set_ylim(-2,1.5)

        ax[-2].set_xlim(-2.1,1.8)
        ax[-2].set_ylim(-2,1.5)
        ax[-1].set_xlim(0,2.1)
        ax[-1].set_ylim(0,2.5)
        plt.savefig(img_dir + f"{direction}_magnetosheath_traversals.svg")

        # Plot histogram of all traversal times
        hist = HistogramPanel(ms_dt, bins = bins, color=color, minmax=True)
        hist.ax_set_params = {
                "title": f"Time between {label} crossings",
                "xlabel": "Time (Hours)",
                "ylabel": "Number of crossings",
                "yscale": "log"
                }
        plots.append(hist)

    dt_hist_bs_mp = plots[0] + plots[1]
    dt_hist_bs_mp.plot(show=False, figsize=(18,16))

    # Print total time spent magnetosheath
    print("Number of days in magnetosheath" , (sum(plots[1]._data) + sum(plots[1]._data))/(24))

plt.show()
