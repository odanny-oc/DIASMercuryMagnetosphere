import os
from pathlib import Path
os.chdir(Path(__file__).resolve().parent)

import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import numpy as np
from astropy.table import QTable, vstack, Column, Table
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
from hermpymod.functions.boundary_models_mod import boundary_fitter, plot_magnetospheric_boundaries
from hermpymod.functions.grazing_angle import get_grazing_angle


home_dir = os.getenv('HOME')
from hermpymod.paths import DATA_DIR
data_dir = DATA_DIR
img_dir = "../../../plots_and_images/"

plt.style.use(img_dir + "presentation.mplstyle")

Zd = c.DIPOLE_OFFSET.to("Mercury Radii")

spice_client = ClientSPICE()

hollman_crossing_list = parse_crossing_list()
hollman_crossing_list["UTC"] = Time(hollman_crossing_list["UTC"]).to_datetime()

hollman_mp_crossing_list = vstack([i for i in hollman_crossing_list if "MP" in i["Label"]])
hollman_bs_crossing_list = vstack([i for i in hollman_crossing_list if "BS" in i["Label"]])


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

    data_sets = [
            (ms_in, r"BS IN \& MP IN", "purple"),
            (ms_out, r"MP OUT \& BS OUT", "pink"),
            ]
    
    plots = []

    id_start =0
    id_end = 50

    for data_type, label, color in data_sets:

        ms_dt = np.array([(i[-1]["UTC"] - i[0]["UTC"]).total_seconds()/3600 for i in data_type])
        
        dt_sorted = np.sort(ms_dt)

        #Take shortest 100
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

        ms_crossings = [vstack(i) for i in ms_short]
        ms_crossings = vstack(ms_crossings)

        # Order ms_crossings in time instead of from shortest traversal to longest
        ms_crossings = ms_crossings[np.argsort(Time(ms_crossings["UTC"]).to_datetime())]

        mp_crossings = ms_crossings[["MP" in i["Label"] for i in ms_crossings]]
        bs_crossings = ms_crossings[["BS" in i["Label"] for i in ms_crossings]]

        if mp_crossings["Label"][0] == "MP_IN":
            mp_indices = np.searchsorted(hollman_mp_crossing_list["UTC"], mp_crossings["UTC"]) + 1
            bs_indices = np.searchsorted(hollman_bs_crossing_list["UTC"], bs_crossings["UTC"]) - 1
        else:
            mp_indices = np.searchsorted(hollman_mp_crossing_list["UTC"], mp_crossings["UTC"]) - 1
            bs_indices = np.searchsorted(hollman_bs_crossing_list["UTC"], bs_crossings["UTC"]) + 1

        indices = np.concatenate([mp_indices, bs_indices])
        indices = sorted(indices)

        adjacent_crossings = hollman_crossing_list[indices]
        with spice_client.KernelPool():

            for idx in shortest_indices:
                ms = data_type[idx]
                time_range = [ms[0]["UTC"].isoformat(),ms[-1]["UTC"].isoformat()]
                plot_all_ephemeris(time_range, color='k', ax=ax, scatter=False, downsampled=False, resolution=10, frame="MSM")

        plot_all_ephemeris([data_type[0][0]["UTC"].isoformat(),data_type[-1][-1]["UTC"].isoformat()], color='none', ax=ax, crossings=ms_crossings, frame="MSM")
        plot_all_ephemeris([data_type[0][0]["UTC"].isoformat(),data_type[-1][-1]["UTC"].isoformat()], color='none', ax=ax, crossings=adjacent_crossings, frame="MSM")

        fig.suptitle(label)

        ms_crossings_abb = ms_crossings.copy()
        adjacent_crossings_abb = adjacent_crossings.copy()
            
        # ms_crossings_abb["X MSO"] = ms_crossings["X MSO"] * np.cos(7 * np.pi/180) + ms_crossings["Y MSO"] * np.sin(7 * np.pi/180)
        # ms_crossings_abb["Y MSO"] = - ms_crossings["X MSO"] * np.sin(7 * np.pi/180) + ms_crossings["Y MSO"] * np.cos(7 * np.pi/180)
        #
        # adjacent_crossings_abb["X MSO"] = adjacent_crossings["X MSO"] * np.cos(7 * np.pi/180) + adjacent_crossings["Y MSO"] * np.sin(7 * np.pi/180)
        # adjacent_crossings_abb["Y MSO"] = - adjacent_crossings["X MSO"] * np.sin(7 * np.pi/180) + adjacent_crossings["Y MSO"] * np.cos(7 * np.pi/180)


        for idx, crossing in enumerate(ms_crossings_abb[:5]):
            crossings_fit = [crossing, adjacent_crossings_abb[idx]]
            if "MP" in crossing["Label"]:
                function = 'Magnetopause'
            else:
                function = 'Bow Shock'

            boundary_parameters = boundary_fitter(crossings_fit, function=function, epsilon=1.02, alpha=0.5)

            crossings_table = vstack(crossings_fit)
            crossings_table = crossings_table[np.argsort(crossings_table["UTC"])]

            time_range =[(Time(crossings_table[0]["UTC"]).to_datetime() - dt.timedelta(hours=1)).isoformat(), (Time(crossings_table[-1]["UTC"]).to_datetime() + dt.timedelta(hours=1)).isoformat()]

            fig, ax = plot_all_ephemeris(time_range, crossings=crossings_table, frame='MSM')

            config = [
            (ax[0], "xy", False),
            (ax[1], "xz", False),
            (ax[2], "yz", False),
            (ax[3], "xy", True),
            ]

            if function == 'Magnetopause':
                for axis, plane, cylindrical in config:
                    plot_magnetospheric_boundaries(axis, plane=plane, frame="MSM", alpha=boundary_parameters[1], sub_solar_magnetopause=boundary_parameters[0], color='red', cylindrical=cylindrical)

            elif function == 'Bow Shock':
                for axis, plane, cylindrical in config:
                    plot_magnetospheric_boundaries(axis, plane=plane, frame="MSM", psi=boundary_parameters[1], p=boundary_parameters[0], color='red', cylindrical=cylindrical)

        
plt.show()
