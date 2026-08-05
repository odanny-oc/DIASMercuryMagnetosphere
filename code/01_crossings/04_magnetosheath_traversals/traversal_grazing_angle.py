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
from hermpymod.functions.ephemeris_downsampler import parse_crossing_list
from hermpymod.functions.grazing_angle import get_grazing_angle


home_dir = os.getenv('HOME')
data_dir = os.path.join(home_dir, '.ephemeris_data/')
img_dir = "../../../plots_and_images/"

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
    dt, ms_out, ms_in = delta_t_magenetosheath(data)

    """
    Plots to make, 
    inbound magnetosheath traversals, coloured by BS grazing angle
    inbound magnetosheath traversals, coloured by MP grazing angle
    outbound magnetosheath traversals, coloured by BS grazing angle
    outbound magnetosheath traversals, coloured by MP grazing angle
    """

    data_sets = [
            (ms_in, r"BS IN \& MP IN", "purple", "BS"),
            (ms_in, r"BS IN \& MP IN", "purple", "MP"),
            (ms_out, r"MP OUT \& BS OUT", "pink", "BS"),
            (ms_out, r"MP OUT \& BS OUT", "pink", "MP"),
            ]
    
    plots = []

    id_start = 0
    id_end = 100

    for data_type, label, color, angle in data_sets:

        # Magnetosheath traversal duration
        ms_dt = np.array([(i[-1]["UTC"] - i[0]["UTC"]).total_seconds()/3600 for i in data_type])
        
        #Take shortest 100 traversals
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

        grazing_angle_vectors = []
        mp_grazing_angle = []
        bs_grazing_angle = []

        # Calculate grazing angle
        with spice_client.KernelPool():
            for crossing in ms_crossings:
                if "MP" in crossing["Label"]:
                    # Returns angle, boundary normal vector, and velocity vector
                    gz_vec = get_grazing_angle(crossing, function="Magnetopause", return_vectors=True)
                    grazing_angle_vectors.append(gz_vec)
                    # Save angle for colour map
                    mp_grazing_angle.append(gz_vec[0])
                if "BS" in crossing["Label"]:
                    gz_vec = get_grazing_angle(crossing, function="Bow Shock", return_vectors=True)
                    grazing_angle_vectors.append(gz_vec)
                    bs_grazing_angle.append(gz_vec[0])

            grazing_angle_vectors = np.array(grazing_angle_vectors, dtype=object)

            # Set colourmap weights
            if angle == "MP":
                values = mp_grazing_angle
            else:
                values = bs_grazing_angle

            cmap = mpl.colormaps['plasma']
            norm = Normalize(vmin=min(values), vmax=max(values))
            

            # Plot shortest 100 traversals, coloured by grazing angle
            for i,idx in enumerate(shortest_indices):
                ms = data_type[idx]
                time_range = [ms[0]["UTC"].isoformat(),ms[-1]["UTC"].isoformat()]
                color=cmap(norm(values[i]))
                plot_all_ephemeris(time_range, color=color, ax=ax, scatter=False, downsampled=False, resolution=10, frame="MSM")


        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, label=angle + ' Grazing angle')

        # Plot magnetic boundaries
        plot_all_ephemeris([data_type[0][0]["UTC"].isoformat(),data_type[-1][-1]["UTC"].isoformat()], color='none', ax=ax, boundaries=True, frame="MSM")

        # Plot normal vectors
        ax[-1].quiver(ms_crossings["X MSO"], np.sqrt(ms_crossings["Y MSO"]**2 + (ms_crossings["Z MSO"] - Zd.value)**2), [arr[0] for arr in grazing_angle_vectors[:,1]], [arr[1] for arr in grazing_angle_vectors[:,1]], color="red")

        # Plot velocity vectors
        ax[-1].quiver(ms_crossings["X MSO"], np.sqrt(ms_crossings["Y MSO"]**2 + (ms_crossings["Z MSO"] - Zd.value)**2), [arr[0] for arr in grazing_angle_vectors[:,2]], [arr[1] for arr in grazing_angle_vectors[:,2]], color='blue')

        fig.suptitle(label)

plt.show()
