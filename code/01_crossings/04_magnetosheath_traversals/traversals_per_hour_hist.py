import os
from pathlib import Path
os.chdir(Path(__file__).resolve().parent)

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import colors
from matplotlib.colors import LogNorm

import numpy as np

import os

from astropy.table import QTable, vstack
from astropy.time import Time
import astropy.units as u

from hermpy.utils.constants import Constants
from hermpy.net.client_spice import ClientSPICE

from hermpymod.classes.panels import PlanarplotPanel, HistogramPanel
from hermpymod.functions.ephemeris_downsampler import parse_crossing_list, abs_r, parse_spice, time_array,  parse_periapsis_data
from hermpymod.functions.downsampled_positional_data import parse_spice_downsampled
from hermpymod.functions.encounters import encounter_finder, parse_encounters_list
from hermpymod.functions.data_per_orbit import orbit_data
from hermpymod.functions.boundary_models_mod import plot_magnetospheric_boundaries
from hermpymod.functions.plot_all import plot_all_ephemeris
from hermpymod.functions.magnetosheath_traversals import magnetosheath_crossings
from hermpymod.paths import DATA_DIR


data_dir = DATA_DIR
img_dir = "../../../plots_and_images/"


Zd = Constants.DIPOLE_OFFSET.to("Mercury Radii")

plt.style.use(img_dir + "presentation.mplstyle")

crossing_data = parse_crossing_list()
crossing_data["UTC"] = Time(crossing_data["UTC"]).to_datetime()

# Could be resolution problem here
ephemeris_data = parse_spice_downsampled()

ephemeris_data["UTC"] = ephemeris_data["UTC"].to_datetime()

residence_time = [(ephemeris_data["UTC"][i+1] - ephemeris_data["UTC"][i]).total_seconds()/3600 for i in range(len(ephemeris_data) - 1)]

# Make residence time same shape as ephemeris_data
residence_time.append(0)

ephemeris_rho = np.sqrt(ephemeris_data["Y MSM"]**2 + ephemeris_data["Z MSM"]**2)

ms_out , ms_in = magnetosheath_crossings(crossing_data)


ms_out_dt = [i[-1]["UTC"] - i[0]["UTC"] for i in ms_out]
ms_in_dt = [i[-1]["UTC"] - i[0]["UTC"] for i in ms_in]

# Take shortest 100 traversals
ms_out = [ms_out[i] for i in np.argsort(ms_out_dt)[:100]]
ms_in = [ms_in[i] for i in np.argsort(ms_in_dt)[:100]]


ms_out_start_times = np.array([i[0]["UTC"] for i in ms_out])
ms_out_endtimes = np.array([i[-1]["UTC"] for i in ms_out])

ms_in_start_times = np.array([i[0]["UTC"] for i in ms_in])
ms_in_endtimes = np.array([i[-1]["UTC"] for i in ms_in])

# Spatial indices of start and end of inward and outward traversals
out_start_idx = np.searchsorted(ephemeris_data["UTC"], ms_out_start_times, side="left")
out_end_idx = np.searchsorted(ephemeris_data["UTC"], ms_out_endtimes, side="right")

in_start_idx = np.searchsorted(ephemeris_data["UTC"], ms_in_start_times, side="left")
in_end_idx = np.searchsorted(ephemeris_data["UTC"], ms_in_endtimes, side="right")


config = [
    (out_start_idx, out_end_idx, "MP_OUT to BS_OUT"),
    (in_start_idx, in_end_idx, "BS_IN to MP_IN"),
]


for start_idx, end_idx, title in config:

    # Bin sizes, 0.01 Mercury Radii ^2 from -7 to 7 RM
    xedges, yedges = np.arange(-7,7, 0.1),np.arange(-7,7, 0.1) 

    # Size of bin matrix (minus one due to NumPy convention)
    nx = len(xedges) - 1
    ny = len(yedges) - 1

    fig = plt.figure()

    # Initialise non-standard sized axis for plotting
    gs = GridSpec(2, 3, figure=fig)

    ax = [
            fig.add_subplot(gs[0, 0]),
            fig.add_subplot(gs[0, 1]),
            fig.add_subplot(gs[0, 2]), 
            fig.add_subplot(gs[1, :]),
            ]

    # Config for each plane and cylindrical coords
    encounter_per_hour_hist_configs = [
            (ephemeris_data["X MSM"], ephemeris_data["Y MSM"], ax[0], "xy", "X MSM", "Y MSM"),
            (ephemeris_data["X MSM"], ephemeris_data["Z MSM"], ax[1], "xz", "X MSM", "Z MSM"),
            (ephemeris_data["Y MSM"], ephemeris_data["Z MSM"], ax[2], "yz", "Y MSM", "Z MSM"),
            (ephemeris_data["X MSM"], ephemeris_rho, ax[3], "cylindrical", "X MSM", r"$\rho$ MSM"),
            ]

    for x, y, axis, label, xname, yname in encounter_per_hour_hist_configs:

        # Change size for Cylindrical coords as they're not symmetric
        if label == "cylindrical":
            yedges = np.arange(0,10, 0.1)
            ny = len(yedges) - 1

        # Initialise bin matrix
        traversal_count = np.zeros((nx,ny))

        # Histogram of 'residence' of MESSENGER during entire mission
        residence_hours, _, _= np.histogram2d(x.value, y.value, bins=[xedges, yedges], weights=residence_time)

        # Find bins that traversals lie in

        for i, (s, e) in enumerate(zip(start_idx, end_idx)):

            if s==e:
                ix = np.digitize(x[s].value, xedges) - 1
                iy = np.digitize(y[s].value, yedges) - 1

            else:
                ix = np.digitize(x[s:e].value, xedges) - 1
                iy = np.digitize(y[s:e].value, yedges) - 1

            # Pair indices as coordinates and remove duplicates (overcounting traversals per bin)
            bin_pairs = np.column_stack([ix, iy])
            unique_bins = np.unique(bin_pairs, axis=0)
            
            # Increase count of bin matrix for positions that traversals crossed through
            traversal_count[unique_bins[:, 0], unique_bins[:, 1]] += 1


        # Find rate of traversals per hour by dividing it (when possible i.e no NaN or 0) by how long MESSENGER was in each bin
        with np.errstate(divide='ignore', invalid='ignore'):
            rate_grid = np.where(residence_hours > 0, traversal_count / residence_hours, np.nan)


        # Plot
        pcm = axis.pcolormesh(xedges, yedges, rate_grid.T, norm=LogNorm(), cmap='viridis')

        if label =="xy":
            fig.colorbar(pcm, ax=ax, label= title + r' Traversals (hours)$^{-1}$ (Normalised by residence)')

        # Plot boundaries and label
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


    fig.suptitle("Traversals per hour in MESSENGER Mission")

plt.show()

