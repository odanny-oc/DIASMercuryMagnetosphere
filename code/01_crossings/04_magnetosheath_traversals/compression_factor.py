import os
import pandas as pd
import numpy as np

from pathlib import Path
os.chdir(Path(__file__).resolve().parent)

import datetime as dt

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import BoundaryNorm
from matplotlib.lines import Line2D

from scipy.stats import binned_statistic_2d

from hermpy.plotting import Panel, TimeseriesPanel
from hermpy.net import ClientSPICE, ClientMESSENGER
import spiceypy as spice
from hermpy.utils import Constants as c

import astropy.units as u
from astropy.time import Time
from astropy.table import QTable, hstack, vstack


from hermpymod.classes.panels import PlanarplotPanel,  HistogramPanel
from hermpymod.functions.plot_all import plot_all_ephemeris
from hermpymod.functions.ephemeris_downsampler import parse_crossing_list
from hermpymod.functions.magnetosheath_traversals import magnetosheath_crossings
from hermpymod.paths import DATA_DIR


data_dir = DATA_DIR
img_dir = "../../../plots_and_images/"

plt.style.use(img_dir + "presentation.mplstyle")

spice_client = ClientSPICE()

Zd = c.DIPOLE_OFFSET.to("Mercury Radii").value

"""
Functions to plot magnetospheric boundaries. Both plot in MSM coordinates. See Winslow et al. paper for details
"""


def shue_surface_msm(theta, phi, Rss, alpha):
    r = shue_model(theta, Rss, alpha)
    X = r * np.cos(theta)
    Y = r * np.sin(theta) * np.cos(phi)
    Z = r * np.sin(theta) * np.sin(phi)
    return X, Y, Z


def conic_surface_msm(theta, phi, p, epsilon, X0=0.5):
    r = conic_section(theta, p, epsilon)
    X = X0 + r * np.cos(theta)
    Y = r * np.sin(theta) * np.cos(phi)
    Z = r * np.sin(theta) * np.sin(phi)
    return X, Y, Z


def shue_model(theta, Rss, alpha):
    return Rss*(2/(1 + np.cos(theta)))**alpha


def conic_section(theta, p, epsilon):
    return (p*epsilon/(1 + epsilon*np.cos(theta)))


hollman_crossing_list = parse_crossing_list()
hollman_crossing_list["UTC"] = Time(hollman_crossing_list["UTC"]).to_datetime()

configs = [
        (hollman_crossing_list, "Hollman")
        ]

for data, name in configs:
    ms_out, ms_in = magnetosheath_crossings(data)

    """
    Plots to make, 
    inbound magnetosheath traversals
    outbound magnetosheath traversals
    """

    data_sets = [
            (ms_in, r"BS IN \& MP IN", "purple", "Inbound"),
            (ms_out, r"MP OUT \& BS OUT", "pink", "Outbound"),
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


        # Stack crossings for plotting
        ms_crossings = [vstack(i) for i in ms_short]
        ms_crossings = vstack(ms_crossings)

        # Split into first and last crossings
        first_crossing = vstack([i[0] for i in data_type])
        second_crossing = vstack([i[-1] for i in data_type])

        # Get radial distance to gauge compression
        radial_distance_first = [np.sqrt(i["X MSO"]**2 + i["Y MSO"]**2 + (i["Z MSO"])**2) for i in first_crossing]
        radial_distance_second = [np.sqrt(i["X MSO"]**2 + i["Y MSO"]**2 + (i["Z MSO"])**2) for i in second_crossing]

        if direction == "Inbound":
            # Project radially the point onto the relevant boundary ( depends on case )

            theta_shue = np.arctan2(np.sqrt(second_crossing["Y MSO"]**2 + (second_crossing["Z MSO"] - Zd)**2), second_crossing["X MSO"])
            phi_shue = np.arctan2((second_crossing["Z MSO"] - Zd), second_crossing["Y MSO"])

            theta_conic = np.arctan2((np.sqrt(first_crossing["Y MSO"]**2 + (first_crossing["Z MSO"] - Zd)**2)), (first_crossing["X MSO"] - 0.5))
            phi_conic = np.arctan2((first_crossing["Z MSO"] - Zd), first_crossing["Y MSO"] )

        else:
            theta_shue = np.arctan2((np.sqrt(first_crossing["Y MSO"]**2 + (first_crossing["Z MSO"] - Zd)**2)), first_crossing["X MSO"])
            phi_shue = np.arctan2((first_crossing["Z MSO"] - Zd), first_crossing["Y MSO"])

            theta_conic = np.arctan2((np.sqrt(second_crossing["Y MSO"]**2 + (second_crossing["Z MSO"] - Zd)**2)),(second_crossing["X MSO"] - 0.5))
            phi_conic = np.arctan2((second_crossing["Z MSO"] - Zd), second_crossing["Y MSO"])

        # Calculate projected boundary position
        mp_pos = shue_surface_msm(theta_shue, phi_shue, Rss=1.45, alpha=0.5)
        bs_pos = conic_surface_msm(theta_conic, phi_conic, p=2.75, epsilon=1.02)

        # Calculate boundary distance to Mercury (MSO coordinates)
        mp_radial_dist = np.sqrt(mp_pos[0]**2 + mp_pos[1]**2 + (mp_pos[2] + Zd)**2)
        bs_radial_dist = np.sqrt(bs_pos[0]**2 + bs_pos[1]**2 + (bs_pos[2] + Zd)**2)
 
        point_size = 0.2
        if direction == "Inbound":

            delta_r_bs = (radial_distance_first - bs_radial_dist)
            delta_r_mp = (radial_distance_second - mp_radial_dist)

        else:

            delta_r_bs = (radial_distance_second - bs_radial_dist)

            delta_r_mp = (radial_distance_first - mp_radial_dist)

        
        if direction == "Inbound":
            hists_config = [
                    (first_crossing, delta_r_bs, bs_radial_dist, "Bow Shock"),
                    (second_crossing, delta_r_mp, mp_radial_dist, "Magnetopause"),
                    ]

        else:
            hists_config = [
                    (first_crossing, delta_r_mp, mp_radial_dist, "Magnetopause"),
                    (second_crossing, delta_r_bs, bs_radial_dist, "Bow Shock"),
                    ]
        

        for crossing, delta_r, boundary, boundary_label in hists_config:
            compression_factor = delta_r / boundary
            df = pd.DataFrame({"UTC" : crossing["UTC"], "Compression": compression_factor}).set_index('UTC')
            df.index = pd.to_datetime(df.index)
            window = "10D"
            x_smooth = df["Compression"].rolling(window).mean()



            fig, ax = plt.subplots()
            ax.plot(x_smooth.index, x_smooth, label=f'Compression factor, Boxcar average ({window})')
            fig.suptitle(f"Compression Factor for {direction} {boundary_label}")
            ax.set_xlabel('Time UTC (YYYY:MM:DD)')
            ax.set_ylabel(f'{boundary_label} Compression')
            ax.axhline(color='k', ls='--')
            ax.legend()


plt.show()
