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
from hermpymod.functions.magnetosheath_traversals import magnetosheath_crossings
from hermpymod.paths import DATA_DIR


data_dir = DATA_DIR
img_dir = "../../../plots_and_images/"

plt.style.use(img_dir + "presentation.mplstyle")

Zd = c.DIPOLE_OFFSET.to("Mercury Radii")

spice_client = ClientSPICE()

hollman_crossing_list = parse_crossing_list()
hollman_crossing_list["UTC"] = Time(hollman_crossing_list["UTC"]).to_datetime()

configs = [
        (hollman_crossing_list, "Hollman")
        ]


for data, name in configs:
    ms_out, ms_in = magnetosheath_crossings(data)

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
        ms_crossings = vstack([vstack(i) for i in ms_short])
        mp_crossings = vstack([i for i in ms_crossings if "MP" in i["Label"]])
        bs_crossings = vstack([i for i in ms_crossings if "BS" in i["Label"]])

        # Calculate grazing angle
        with spice_client.KernelPool():
            mp_rho = np.sqrt(mp_crossings["Y MSO"]**2 + (mp_crossings["Z MSO"] - Zd.value)**2)
            mp_pos = np.array([mp_crossings["X MSO"], mp_rho])

            # Returns angle, boundary normal vector, and velocity vector
            mp_grazing_angle = get_grazing_angle(time = mp_crossings["UTC"], position=mp_pos, function="Magnetopause", return_vectors=True)
            # Save angle for colour map

            bs_rho = np.sqrt(bs_crossings["Y MSO"]**2 + (bs_crossings["Z MSO"] - Zd.value)**2)
            bs_pos = np.array([bs_crossings["X MSO"], bs_rho])
            bs_grazing_angle = get_grazing_angle(time=bs_crossings["UTC"], position=bs_pos, function="Bow Shock", return_vectors=True)


            print(bs_grazing_angle)

            # Set colourmap weights
            if angle == "MP":
                values = mp_grazing_angle[0]
            else:
                values = bs_grazing_angle[0]

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
        ax[-1].quiver(mp_crossings["X MSO"], mp_rho, mp_grazing_angle[1][0], mp_grazing_angle[1][1], color="blue")

        ax[-1].quiver(bs_crossings["X MSO"], bs_rho, bs_grazing_angle[1][0], bs_grazing_angle[1][1], color="blue")

        # Plot velocity vectors
        ax[-1].quiver(mp_crossings["X MSO"], mp_rho, mp_grazing_angle[2][0], mp_grazing_angle[2][1], color="red")

        ax[-1].quiver(bs_crossings["X MSO"], bs_rho, bs_grazing_angle[2][0], bs_grazing_angle[2][1], color="red")

        fig.suptitle(label)

plt.show()
