import astropy.units as u
import matplotlib.pyplot as plt

from matplotlib.gridspec import GridSpec
from matplotlib.colors import LogNorm


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

from hermpy.utils.constants import Constants

import os

from hermpymod.classes.panels import PlanarplotPanel,  HistogramPanel
from hermpymod.functions.plot_all import plot_all_ephemeris
from hermpymod.functions.ephemeris_downsampler import parse_crossing_list
from hermpymod.functions.downsampled_positional_data import parse_spice_downsampled
from hermpymod.functions.boundary_models_mod import plot_magnetospheric_boundaries


home_dir = os.getenv('HOME')
data_dir = os.path.join(home_dir, '.ephemeris_data/')
img_dir = "../../plots_and_images/"

plt.style.use(img_dir + "presentation.mplstyle")

Zd = Constants.DIPOLE_OFFSET.to("Mercury Radii")

hollman_crossing_list = parse_crossing_list()

# MSM coords for crossings
hollman_crossing_list["X MSM"] = hollman_crossing_list["X MSO"]
hollman_crossing_list["Y MSM"] = hollman_crossing_list["Y MSO"]
hollman_crossing_list["Z MSM"] = hollman_crossing_list["Z MSO"] + Zd

hollman_crossing_list["UTC"] = Time(hollman_crossing_list["UTC"]).to_datetime()

ephemeris_data = parse_spice_downsampled()
# This is crucial for vectorisation of np.searchsorted later
ephemeris_data["UTC"] = ephemeris_data["UTC"].to_datetime()

ephemeris_rho = np.sqrt(ephemeris_data["Y MSM"]**2 + ephemeris_data["Z MSM"]**2)


"""
Function that finds all crossings of BS_IN -> MP_IN and MP_OUT -> BS_OUT
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

    data_sets = [
            (ms_in, r"BS IN -\> MP IN", "purple"),
            (ms_out, r"MP OUT -\> BS OUT", "pink"),
            ]
    

    for data_type, label, color in data_sets:

        # Find length of times that MESSENGER was in the magnetosheath
        ms_dt = np.array([(i[-1]["UTC"] - i[0]["UTC"]).total_seconds()/3600 for i in data_type])

        # Take shortest examples (everything less than the average)
        short_ms_mask = ms_dt < np.average(ms_dt)
        ms_dt_short = ms_dt[short_ms_mask]

        ms_short = []

        for idx, mask in enumerate(short_ms_mask):
            if mask:
                ms_short.append(data_type[idx])

        ms_short_start_times = np.array([i[0]["UTC"] for i in ms_short])
        ms_short_end_times = np.array([i[-1]["UTC"] for i in ms_short])

        # Find indices in positional data where these trajectories occur
        start_idx = np.searchsorted(ephemeris_data["UTC"], ms_short_start_times, side="left")
        end_idx = np.searchsorted(ephemeris_data["UTC"], ms_short_end_times, side="right")

        # Print how many short magnetosheath traversals there are
        print("Number of 'short' magnetosheath traversals", len(ms_short))
        
        # Define bin sizes
        xedges, yedges = np.arange(-7,7, 0.1),np.arange(-7,7, 0.1) 

        nx = len(xedges) - 1
        ny = len(yedges) - 1

        fig = plt.figure()
        gs = GridSpec(2, 3, figure=fig)

        ax = [
        fig.add_subplot(gs[0, 0]),
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[0, 2]), 
        fig.add_subplot(gs[1, :]),
        ]
        
        magnetosheath_configs = [
            (ephemeris_data["X MSM"], ephemeris_data["Y MSM"], ax[0], "xy", "X MSM", "Y MSM"),
            (ephemeris_data["X MSM"], ephemeris_data["Z MSM"], ax[1], "xz", "X MSM", "Z MSM"),
            (ephemeris_data["Y MSM"], ephemeris_data["Z MSM"], ax[2], "yz", "Y MSM", "Z MSM"),
            (ephemeris_data["X MSM"], ephemeris_rho, ax[3], "cylindrical", "X MSM", r"$\rho$ MSM"),
            ]

        fig.suptitle(label)


        for x, y, axis, label, xname, yname in magnetosheath_configs:

            # Change bin for cylindrical plot (cylindrical always non-negative)
            if label == "cylindrical":
                yedges = np.arange(0,10, 0.1)
                ny = len(yedges) - 1

            dt_of_trajectory = np.empty((nx,ny), dtype=object)

            # Create bin matrix that is an array of empty lists
            for i in range(nx):
                for j in range(ny):
                    dt_of_trajectory[i,j] = []

            for idx, (s, e) in enumerate(zip(start_idx, end_idx)):

                if s==e:
                    ix = np.digitize(x[s].value, xedges) - 1
                    iy = np.digitize(y[s].value, yedges) - 1

                else:
                    ix = np.digitize(x[s:e].value, xedges) - 1
                    iy = np.digitize(y[s:e].value, yedges) - 1


                # Find coordinates of bins that the traversal passed through
                bin_pairs = np.column_stack([ix, iy])
                unique_bins = np.unique(bin_pairs, axis=0)
                
                # Add traversal duration to bins in which were passed through
                for bin in dt_of_trajectory[unique_bins[:, 0], unique_bins[:, 1]]:
                    bin.append(ms_dt_short[idx])

            
            # For non-empty bins, take bin weighting to be the average duration of the magnetosheath traversal
            for i in range(nx):
                for j in range(ny):
                    if dt_of_trajectory[i,j] == []:
                        dt_of_trajectory[i,j] = np.nan
                    else:
                        dt_of_trajectory[i,j] = np.average(dt_of_trajectory[i,j])

            # Change type to float for plotting
            dt_of_trajectory = dt_of_trajectory.astype(float)

            # Plot with reversed colormap (interested in shortest traversals)
            pcm = axis.pcolormesh(xedges, yedges, dt_of_trajectory.T, norm=LogNorm(), cmap='viridis_r')

            # Plot magnetic boundaries
            if label == "cylindrical":
                plot_magnetospheric_boundaries(axis, frame='MSM', sub_solar_magnetopause=1.4, cylindrical=True, add_legend=True, zorder=5)
                axis.set_xlim(-8,6)
                axis.set_ylim(0,10)
                axis.legend()
            else:
                plot_magnetospheric_boundaries(axis, frame='MSM', sub_solar_magnetopause=1.4, plane=label, add_legend=False, zorder=5)
                axis.set_xlim(-8,8)
                axis.set_ylim(-8,8)

            axis.set_xlabel(f"{xname} Mercury Radii")
            axis.set_ylabel(f"{yname} Mercury Radii")

            if label =="xy":
                fig.colorbar(pcm, ax=ax, label=f'Average traversal duration of magnetosheath (hours)')

plt.show()
