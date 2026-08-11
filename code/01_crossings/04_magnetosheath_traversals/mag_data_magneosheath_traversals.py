import os
from pathlib import Path
os.chdir(Path(__file__).resolve().parent)

import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import matplotlib as mpl

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

from hermpymod.classes.panels import PlanarplotPanel,  HistogramPanel, plot_magnetospheric_boundaries
from hermpymod.functions.plot_all import plot_all_ephemeris
from hermpymod.functions.mag_data_plotter import mag_data_plotter
from hermpymod.functions.ephemeris_downsampler import parse_crossing_list
from hermpymod.functions.grazing_angle import get_grazing_angle


home_dir = os.getenv('HOME')
from hermpymod.paths import DATA_DIR
data_dir = DATA_DIR
img_dir = "../../../plots_and_images/"

mpl.use('TkAgg')

plt.style.use(img_dir + "presentation.mplstyle")

Zd = c.DIPOLE_OFFSET.to("Mercury Radii")

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


for data, name in configs:
    _, ms_out, ms_in = delta_t_magenetosheath(data)

    """
    Plots to make, 
    inbound magnetosheath traversals, coloured by BS grazing angle
    outbound magnetosheath traversals, coloured by MP grazing angle
    """

    data_sets = [
            (ms_in, r"BS IN \& MP IN", "purple", "BS", "inbound"),
            (ms_out, r"MP OUT \& BS OUT", "pink", "MP", "outbound"),
            ]
    
    plots = []

    id_start = 0
    id_end = 2

    for data_type, label, color, angle, direction in data_sets:

        # Magnetosheath traversal duration
        ms_dt = np.array([(i[-1]["UTC"] - i[0]["UTC"]).total_seconds()/3600 for i in data_type])
        
        #Take shortest 100 traversals
        shortest_indices = np.argsort(ms_dt)[id_start:id_end]

        ms_short = []

        for idx in shortest_indices:
            ms_short.append(data_type[idx])


        # Stack crossings for plotting
        ms_crossings = [vstack(i) for i in ms_short]
        ms_crossings = vstack(ms_crossings)

        grazing_angle_vectors = []
        mp_grazing_angle = []
        bs_grazing_angle = []

        tw = dt.timedelta(hours=3)


        for idx, ms in enumerate(ms_short):
            fig, ax = mag_data_plotter([ms[0]["UTC"] - tw, ms[-1]["UTC"] + tw], crossings=True, label=f" Traversal Duration {(ms[-1]["UTC"] - ms[0]["UTC"]).total_seconds():.2f}s", zoom=[ms[0]["UTC"] - dt.timedelta(minutes=5), ms[-1]["UTC"] + dt.timedelta(minutes=5)])
            # fig.suptitle(label)

                # Fix fontsize for subfigure
            for axis in ax:
                axis.set_title(axis.get_title(), fontsize=18)
                axis.xaxis.label.set_size(14)
                axis.yaxis.label.set_size(14)
                for label in axis.get_xticklabels() + axis.get_yticklabels():
                    label.set_fontsize(12)

                legend = axis.get_legend()
                for text in legend.get_texts():
                    text.set_fontsize(14)

            plt.savefig(img_dir + f"ms_traversal_mag_data_{idx}_{direction}.svg")
            
            plt.show()
