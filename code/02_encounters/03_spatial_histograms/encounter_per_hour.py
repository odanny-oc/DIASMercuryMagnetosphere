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


home_dir = os.getenv('HOME')
data_dir = os.path.join(home_dir, '.ephemeris_data/')
img_dir = "../../../plots_and_images/"


Zd = Constants.DIPOLE_OFFSET.to("Mercury Radii")

plt.style.use(img_dir + "presentation.mplstyle")

crossing_data = parse_crossing_list()

ephemeris_data = parse_spice_downsampled()

ephemeris_data["UTC"] = ephemeris_data["UTC"].to_datetime()

residence_time = [(ephemeris_data["UTC"][i+1] - ephemeris_data["UTC"][i]).total_seconds()/3600 for i in range(len(ephemeris_data) - 1)]

# Make residence time same shape as ephemeris_data
residence_time.append(0)

ephemeris_rho = np.sqrt(ephemeris_data["Y MSM"]**2 + ephemeris_data["Z MSM"]**2)

hollman_encounters_data = parse_encounters_list()

encounter_bs = vstack([i for i in hollman_encounters_data if "BS" in i["Label"]])
encounter_mp = vstack([i for i in hollman_encounters_data if "MP" in i["Label"]])

bs_encounter_start_time = Time(encounter_bs["Time Start"]).to_datetime()
bs_encounter_end_time = Time(encounter_bs["Time End"]).to_datetime()

mp_encounter_start_time = Time(encounter_mp["Time Start"]).to_datetime()
mp_encounter_end_time = Time(encounter_mp["Time End"]).to_datetime()

# Spatial indices of start and end of BS and MP encounters
bs_start_idx = np.searchsorted(ephemeris_data["UTC"], bs_encounter_start_time, side="left")
bs_end_idx = np.searchsorted(ephemeris_data["UTC"], bs_encounter_end_time, side="right")

mp_start_idx = np.searchsorted(ephemeris_data["UTC"], mp_encounter_start_time, side="left")
mp_end_idx = np.searchsorted(ephemeris_data["UTC"], mp_encounter_end_time, side="right")


config = [
    (bs_start_idx, bs_end_idx, "Bow Shock"),
    (mp_start_idx, mp_end_idx, "Magnetopause"),
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
        encounter_count = np.zeros((nx,ny))

        # Histogram of 'residence' of MESSENGER during entire mission
        residence_hours, _, _= np.histogram2d(x.value, y.value, bins=[xedges, yedges], weights=residence_time)

        hist_fig, hist_ax = plt.subplots()

        # Plot residence histogram (transpose as first index (x) is a row not column)
        hist_2d = hist_ax.pcolormesh(xedges, yedges, residence_hours.T, norm=LogNorm(), cmap='viridis')
        hist_ax.set_xlabel(xname)
        hist_ax.set_ylabel(yname)

        hist_fig.colorbar(hist_2d, ax=hist_ax, label=r'MESSENGER Residence [hours]')

        
        # Find bins that encounters lie in

        for i, (s, e) in enumerate(zip(start_idx, end_idx)):

            # Consider encounters of single crossing
            if s==e:
                ix = np.digitize(x[s].value, xedges) - 1
                iy = np.digitize(y[s].value, yedges) - 1

            else:
                ix = np.digitize(x[s:e].value, xedges) - 1
                iy = np.digitize(y[s:e].value, yedges) - 1

            # Pair indices as coordinates and remove duplicates (overcounting encounters per bin)
            bin_pairs = np.column_stack([ix, iy])
            unique_bins = np.unique(bin_pairs, axis=0)
            
            # Increase count of bin matrix for positions that encounter crossed through
            encounter_count[unique_bins[:, 0], unique_bins[:, 1]] += 1


        # Find rate of encounters per hour by dividing it (when possible i.e no NaN or 0) by how long MESSENGER was in each bin
        with np.errstate(divide='ignore', invalid='ignore'):
            rate_grid = np.where(residence_hours > 0, encounter_count / residence_hours, np.nan)


        # Plot
        pcm = axis.pcolormesh(xedges, yedges, rate_grid.T, norm=LogNorm(), cmap='viridis')

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


        if label =="xy":
            fig.colorbar(pcm, ax=ax, label= title + r' Encounters (hours)$^{-1}$ (Normalised by residence)')

    fig.suptitle("Encounters per hour in MESSENGER Mission")

plt.show()
