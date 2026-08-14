import os
import pandas as pd
from pathlib import Path
os.chdir(Path(__file__).resolve().parent)

import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import BoundaryNorm

from scipy.stats import binned_statistic_2d

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
 
        fig, ax = plt.subplots(2)
        point_size = 0.2
        if direction == "Inbound":
            ax[0].scatter(first_crossing["UTC"],radial_distance_first, label="Last BS IN", s=point_size)
            ax[0].scatter(first_crossing["UTC"],bs_radial_dist, label="Average BS boundary distance from Mercury, in radial line with MESSENGER", s=point_size)
            delta_r_bs = (radial_distance_first - bs_radial_dist)
            ax[0].scatter(first_crossing["UTC"], delta_r_bs, label="Distance from MESSENGER to average BS boundary", s=point_size)

            ax[1].scatter(second_crossing["UTC"],radial_distance_second,label= "First MP IN", s=point_size)
            ax[1].scatter(second_crossing["UTC"],mp_radial_dist, label="Average MP boundary distance from Mercury, in radial line with MESSENGER" , s=point_size)
            delta_r_mp = (radial_distance_second - mp_radial_dist)
            ax[1].scatter(second_crossing["UTC"], delta_r_mp, label ="Distance from MESSENGER to average MP boundary", s=point_size)
        else:
            ax[0].scatter(second_crossing["UTC"],radial_distance_second, label = "First BS OUT", s=point_size)
            ax[0].scatter(second_crossing["UTC"],bs_radial_dist, label="Average BS distance from Mercury, in radial line with MESSENGER", s=point_size)

            delta_r_bs = (radial_distance_second - bs_radial_dist)
            ax[0].scatter(second_crossing["UTC"], delta_r_bs, label="Distance from average BS boundary to MESSENGER", s=point_size)

            ax[1].scatter(first_crossing["UTC"],radial_distance_first, label="Last MP OUT", s=point_size)
            ax[1].scatter(first_crossing["UTC"],mp_radial_dist, label="Average MP distance from Mercury, in radial line with MESSENGER", s=point_size)

            delta_r_mp = (radial_distance_first - mp_radial_dist)
            ax[1].scatter(first_crossing["UTC"], delta_r_mp, label="Distance from average MP boundary to MESSENGER", s=point_size)

        for axis in ax:
            axis.axhline(color='k', ls='--')

        fig.suptitle(direction + " traversals of magnetosheath, distances from Mercury")
        fig.supylabel("Distance (Mercury Radii)")
        fig.supxlabel("Time UTC (YYYY:MM:DD)")

        ax[0].legend()
        ax[1].legend()
        
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
            df = pd.DataFrame({"UTC" : crossing["UTC"], "Distance": delta_r}).set_index('UTC')
            df.index = pd.to_datetime(df.index)
            window = "10D"
            x_smooth = df["Distance"].rolling(window).mean()

            fig, ax = plt.subplots()
            ax.plot(x_smooth.index, x_smooth, label=f'Compression distance Boxcar averaged ({window})')
            fig.suptitle(f"Compression distance for {direction} {boundary_label}")
            ax.set_xlabel('Time UTC (YYYY:MM:DD)')
            ax.set_ylabel('Distance (Mercury Radii)')
            ax.axhline(color='k', ls='--')
            ax.legend()

            compression_factor = delta_r / boundary

            theta_bins = np.linspace(0, np.pi, 20)
            compression_bins = np.linspace(-1, 1, 100)
            

            bin_indices = np.digitize(theta_shue, theta_bins)
            stat, t_edges, theta_edges, binnumber = binned_statistic_2d(
            compression_factor,theta_shue,ms_dt,
            statistic='mean',
            bins=[compression_bins, theta_bins]
        )

            
            bounds = np.concatenate([np.linspace(0, 1, 100), np.linspace(1, max(ms_dt), 5)])
            bounds = np.unique(bounds)
            print(bounds)
            n_bins = len(bounds) - 1

            cmap = plt.get_cmap('viridis_r', n_bins)
            norm = BoundaryNorm(bounds, cmap.N)


            fig, ax = plt.subplots(figsize=(10, 5))
            pcm = ax.pcolormesh(
                t_edges, theta_edges, stat.T,
                cmap=cmap,
                norm=norm,
                shading='auto'
            )
            fig.colorbar(pcm, ax=ax, label='Mean Traversal Duration (hours)', boundaries=bounds)
            fig.suptitle(f"{direction} {boundary_label} Sheath Traversal duration, positional dependence")

            ax.set_xlabel('Compression Factor')
            ax.set_ylabel(r'$\theta$ Radians')


plt.show()
