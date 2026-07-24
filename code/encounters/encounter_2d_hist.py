import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import colors

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

home_dir = os.getenv('HOME')
data_dir = os.path.join(home_dir, '.ephemeris_data/')
img_dir = "../../plots_and_images/"


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

plt.style.use("../../plots_and_images/presentation.mplstyle")

crossing_data = parse_crossing_list()

bins = np.arange(0,7, 1)

hollman_encounters_data = encounter_finder(crossing_data)


average_encounter = []


for encounter in hollman_encounters_data:

    encounter["X MSM"] = encounter["X MSO"]
    encounter["Y MSM"] = encounter["Y MSO"]
    encounter["Z MSM"] = encounter["Z MSO"] + Zd
    
    if len(encounter) == 1:
        x_pos, y_pos, z_pos = encounter["X MSM"][0], encounter["Y MSM"][0], encounter["Z MSM"][0]
        time = encounter["UTC"][0]

    else:
        _, midpoint, time, _ = arclength_midpoint(encounter["UTC"] ,encounter["X MSM"],encounter["Y MSM"] ,encounter["Z MSM"])

        x_pos, y_pos, z_pos = midpoint

    label = encounter["Label"][0][:2] + encounter["Trajectory Direction"][0][0]
    average_encounter.append([time,x_pos,y_pos,z_pos, label])


average_encounter = np.array(average_encounter)

positions =[average_encounter[:,1].astype(float), average_encounter[:,2].astype(float), average_encounter[:,3].astype(float)] 

distance = abs_r(positions)

crossings_per_encounter = [len(i) for i in hollman_encounters_data]
crossings_per_encounter = np.array(crossings_per_encounter)

average_encounter = QTable({
    "UTC": average_encounter[:,0],
    "|R|": distance,
    "X MSM": positions[0],
    "Y MSM": positions[1],
    "Z MSM": positions[2],
    "Rho MSM": np.sqrt(positions[1]**2 + positions[2]**2),
    "Label": average_encounter[:,4],
    "Number of Crossing": crossings_per_encounter,
    })


bs_encounters = vstack([i for i in average_encounter if "BS" in i["Label"]])
mp_encounters = vstack([i for i in average_encounter if "MP" in i["Label"]])

bs_mask = ["BS" in item for item in average_encounter["Label"]]
mp_mask = ["MP" in item for item in average_encounter["Label"]]

type_config = [
        (bs_encounters, "Bow Shock", bs_mask),
        (mp_encounters, "Magnetopause", mp_mask),
        ]


for encounter, title, mask in type_config:

    fig = plt.figure()

    gs = GridSpec(2, 3, figure=fig)

    ax = [
            fig.add_subplot(gs[0, 0]),
            fig.add_subplot(gs[0, 1]),
            fig.add_subplot(gs[0, 2]), 
            fig.add_subplot(gs[1, :]),
            ]

    plot_configs = [
            ("xy", encounter["X MSM"], encounter["Y MSM"], ax[0]),
            ("xz", encounter["X MSM"], encounter["Z MSM"], ax[1]),
            ("yz", encounter["Y MSM"], encounter["Z MSM"], ax[2]),
            ("Cylindrical", encounter["X MSM"], encounter["Rho MSM"], ax[3]),
            ]

    crossings_per_type = crossings_per_encounter[mask]

    for label, xaxis, yaxis, axis in plot_configs:

        sc = axis.hist2d(xaxis, yaxis, weights=crossings_per_type, bins=50, norm=colors.LogNorm())

        
        if label == "Cylindrical":
            plot_magnetospheric_boundaries(axis, frame='MSM', sub_solar_magnetopause=1.4, cylindrical=True, add_legend=True, zorder=5)
            axis.set_xlim(-8,6)
            axis.set_ylim(0,10)
            axis.legend()
        else:
            plot_magnetospheric_boundaries(axis, frame='MSM', sub_solar_magnetopause=1.4, plane=label, add_legend=False, zorder=5)
            axis.set_xlim(-8,8)
            axis.set_ylim(-8,8)

        axis.set_xlabel(f"{xaxis.name} Mercury Radii")
        axis.set_ylabel(f"{yaxis.name} Mercury Radii")


        if label =="xy":
            fig.colorbar(sc[-1], ax=ax, label=r'Number of crossings')


    fig.suptitle(title + " Encounters in MESSENGER Mission")

plt.show()
