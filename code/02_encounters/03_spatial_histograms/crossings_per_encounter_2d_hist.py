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
from hermpymod.functions.ephemeris_downsampler import parse_crossing_list, abs_r, parse_spice, time_array
from hermpymod.functions.encounters import encounter_finder, parse_encounters_list
from hermpymod.functions.data_per_orbit import orbit_data
from hermpymod.functions.ephemeris_downsampler import parse_periapsis_data
from hermpymod.functions.boundary_models_mod import plot_magnetospheric_boundaries
from hermpymod.functions.plot_all import plot_all_ephemeris
from hermpymod.functions.downsampled_positional_data import parse_spice_downsampled

home_dir = os.getenv('HOME')
from hermpymod.paths import DATA_DIR
data_dir = DATA_DIR
img_dir = "../../../plots_and_images/"


def arclength_midpoint(t, x, y, z=None):
    """
    Find the midpoint (by arc length) of a curve given arrays of coordinates.
    x, y, (z) can be in any consistent units (e.g. km, or Mercury radii).
    Returns the index of the point closest to the midpoint, the interpolated
    midpoint coordinates, and the cumulative arc length array.
    """
    coords = np.column_stack([x, y] if z is None else [x, y, z])
    
    # distance between consecutive points
    diffs = np.diff(coords, axis=0)
    seg_lengths = np.linalg.norm(diffs, axis=1)
    
    # cumulative arc length, starting at 0
    cumlen = np.concatenate([[0], np.cumsum(seg_lengths)])
    total_length = cumlen[-1]
    half_length = total_length / 2
    
    # index of last point before the midpoint
    idx = np.searchsorted(cumlen, half_length) - 1
    idx = np.clip(idx, 0, len(coords) - 2)
    
    # linear interpolation between coords[idx] and coords[idx+1]
    seg_frac = (half_length - cumlen[idx]) / seg_lengths[idx]
    midpoint = coords[idx] + seg_frac * diffs[idx]
    time = t[idx] + (t[idx+1] - t[idx])/2
    
    return idx, midpoint, time, cumlen


Zd = Constants.DIPOLE_OFFSET.to("Mercury Radii")

plt.style.use(img_dir + "presentation.mplstyle")

crossing_data = parse_crossing_list()

bins = np.arange(0,7, 1)

hollman_encounters_data = encounter_finder(crossing_data)
ephemeris_data = parse_spice_downsampled()
ephemeris_data["UTC"] = ephemeris_data["UTC"].to_datetime()

ephemeris_data["Rho MSM"] = np.sqrt(ephemeris_data["Y MSM"]**2 + ephemeris_data["Z MSM"]**2)

average_encounter = []

encounter_list = parse_encounters_list()
encounter_times_start = Time(encounter_list["Time Start"]).to_datetime()
encounter_times_end = Time(encounter_list["Time End"]).to_datetime()

encounter_index_start = np.searchsorted(ephemeris_data["UTC"], encounter_times_start)
encounter_index_end = np.searchsorted(ephemeris_data["UTC"], encounter_times_end, "right")

encounter_indices = list(zip(encounter_index_start.tolist(), encounter_index_end.tolist()))

encounter_indices = np.array(encounter_indices)

crossings_per_encounter = [len(i) for i in hollman_encounters_data]
crossings_per_encounter = np.array(crossings_per_encounter)

bs_mask = ["BS" in item for item in encounter_list["Label"]]
mp_mask = ["MP" in item for item in encounter_list["Label"]]

type_config = [
        ("Bow Shock", bs_mask),
        ("Magnetopause", mp_mask),
        ]


for title, mask in type_config:

    start_idx, end_idx = encounter_indices[mask].T

    fig = plt.figure()

    gs = GridSpec(2, 3, figure=fig)

    ax = [
            fig.add_subplot(gs[0, 0]),
            fig.add_subplot(gs[0, 1]),
            fig.add_subplot(gs[0, 2]), 
            fig.add_subplot(gs[1, :]),
            ]


    plot_configs = [
            ("xy", ephemeris_data["X MSM"], ephemeris_data["Y MSM"], ax[0], "X MSM", "Y MSM"),
            ("xz", ephemeris_data["X MSM"], ephemeris_data["Z MSM"], ax[1], "X MSM", "Z MSM"),
            ("yz", ephemeris_data["Y MSM"], ephemeris_data["Z MSM"], ax[2], "Y MSM", "Z MSM"),
            ("Cylindrical", ephemeris_data["X MSM"], ephemeris_data["Rho MSM"], ax[3], "X MSM", r"$\rho$ MSM"),
            ]

    crossings_per_type = crossings_per_encounter[mask]

    for label, xaxis, yaxis, axis, xlab, ylab in plot_configs:


        if label != "Cylindrical":
            xedges, yedges = np.arange(-7, 7, 0.1), np.arange(-7, 7, 0.1)
        else:
            xedges, yedges = np.arange(-7, 7, 0.1), np.arange(0, 10, 0.1)

        nx, ny = len(xedges) - 1, len(yedges) - 1

        
        number_of_crossings = np.empty((nx,ny), dtype=object)

        # Create bin matrix that is an array of empty lists
        for i in range(nx):
            for j in range(ny):
                number_of_crossings[i,j] = []

        for idx, (s, e) in enumerate(zip(start_idx, end_idx)):

            if s==e:
                ix = np.digitize(xaxis[s].value, xedges) - 1
                iy = np.digitize(yaxis[s].value, yedges) - 1

            else:
                ix = np.digitize(xaxis[s:e].value, xedges) - 1
                iy = np.digitize(yaxis[s:e].value, yedges) - 1


            # Find coordinates of bins that the traversal passed through
            bin_pairs = np.column_stack([ix, iy])
            unique_bins = np.unique(bin_pairs, axis=0)
            
            # Add traversal duration to bins in which were passed through
            for bin in number_of_crossings[unique_bins[:, 0], unique_bins[:, 1]]:
                bin.append(crossings_per_type[idx])

        
        # For non-empty bins, take bin weighting to be the average duration of the magnetosheath traversal
        for i in range(nx):
            for j in range(ny):
                if number_of_crossings[i,j] == []:
                    number_of_crossings[i,j] = np.nan
                else:
                    number_of_crossings[i,j] = np.average(number_of_crossings[i,j])


        number_of_crossings = number_of_crossings.astype(float)

        sc = axis.pcolormesh(xedges, yedges, number_of_crossings.T, norm=LogNorm(), cmap='viridis')

        
        if label == "Cylindrical":
            plot_magnetospheric_boundaries(axis, frame='MSM', sub_solar_magnetopause=1.4, cylindrical=True, add_legend=True, zorder=5)
            axis.set_xlim(-8,6)
            axis.set_ylim(0,10)
            axis.legend()
        else:
            plot_magnetospheric_boundaries(axis, frame='MSM', sub_solar_magnetopause=1.4, plane=label, add_legend=False, zorder=5)
            axis.set_xlim(-8,8)
            axis.set_ylim(-8,8)

        axis.set_xlabel(xlab + " Mercury Radii")
        axis.set_ylabel(ylab + " Mercury Radii")


        if label =="xy":
            fig.colorbar(sc, ax=ax, label=r'Number of crossings')


    fig.suptitle(title + " Encounters in MESSENGER Mission")

plt.show()
