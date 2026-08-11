import os
from pathlib import Path
os.chdir(Path(__file__).resolve().parent)

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import BoundaryNorm

import numpy as np

import os

from astropy.table import QTable, vstack
from astropy.time import Time
import astropy.units as u
from hermpy.utils.constants import Constants

from hermpymod.classes.panels import PlanarplotPanel, HistogramPanel
from hermpymod.functions.ephemeris_downsampler import parse_crossing_list, abs_r
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

encounter_list = parse_encounters_list()
encounter_times_start = Time(encounter_list["Time Start"]).to_datetime()
encounter_times_end = Time(encounter_list["Time End"]).to_datetime()

encounter_index_start = np.searchsorted(ephemeris_data["UTC"], encounter_times_start)
encounter_index_end = np.searchsorted(ephemeris_data["UTC"], encounter_times_end, "right")

encounter_indices = list(zip(encounter_index_start.tolist(), encounter_index_end.tolist()))

average_encounter = []


for idx, (s,e) in enumerate(encounter_indices):

    if len(ephemeris_data[s:e]) == 1 or len(ephemeris_data[s:e]) == 0:
        # Convert to MSM coordinates
        x_pos, y_pos, z_pos = ephemeris_data["X MSM"][s].value, ephemeris_data["Y MSM"][s].value, ephemeris_data["Z MSM"][s].value
        time = ephemeris_data["UTC"][s]

    else:
        _, midpoint, time, _ = arclength_midpoint(ephemeris_data["UTC"][s:e].value ,ephemeris_data["X MSM"][s:e].value,ephemeris_data["Y MSM"][s:e].value ,ephemeris_data["Z MSM"][s:e].value)

        x_pos, y_pos, z_pos = midpoint

    label = encounter_list["Label"][idx]
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


encounter_err = []


for encounter in hollman_encounters_data:
    x_err = abs(encounter[0]["X MSO"] - encounter[-1]["X MSO"])
    y_err = abs(encounter[0]["Y MSO"] - encounter[-1]["Y MSO"])
    z_err = abs(encounter[0]["Z MSO"] - encounter[-1]["Z MSO"])
    rho_err = abs(np.sqrt(encounter[0]["Y MSO"]**2 + (encounter[0]["Z MSO"] + Zd.value)**2) - np.sqrt(encounter[-1]["Y MSO"]**2 + (encounter[-1]["Z MSO"] + Zd.value)**2))
    
    encounter_err.append([x_err, y_err, z_err, rho_err])


encounter_err = np.array(encounter_err)

bs_encounters = ["BS" in i for i in average_encounter["Label"]]
mp_encounters = ["MP" in i for i in average_encounter["Label"]]

average_encounter_type = [
        (bs_encounters, "Bow Shock"), 
        (mp_encounters, "Magnetopause"),
                                       ]
    

for mask, title in average_encounter_type:

    fig = plt.figure()

    gs = GridSpec(2, 3, figure=fig)

    ax = [
            fig.add_subplot(gs[0, 0]),
            fig.add_subplot(gs[0, 1]),
            fig.add_subplot(gs[0, 2]), 
            fig.add_subplot(gs[1, :]),
            ]

    encounters = average_encounter[mask]
    encoun_err = encounter_err[mask]

    plot_configs = [
            ("xy", encounters["X MSM"], encounters["Y MSM"], encoun_err[:,0], encoun_err[:,1], ax[0]),
            ("xz", encounters["X MSM"], encounters["Z MSM"], encoun_err[:,0], encoun_err[:,2], ax[1]),
            ("yz", encounters["Y MSM"], encounters["Z MSM"], encoun_err[:,1], encoun_err[:,2], ax[2]),
            ("Cylindrical", encounters["X MSM"], encounters["Rho MSM"], encoun_err[:,0], encoun_err[:,3], ax[3]),
            ]

    crossings_per_type = crossings_per_encounter[mask]

    bounds = np.concatenate([np.linspace(1, 10, 10), np.linspace(11, max(crossings_per_type), 5)]).astype(int)
    bounds = np.unique(bounds)
    print(bounds)
    n_bins = len(bounds) - 1

    cmap = plt.get_cmap('viridis', n_bins)
    norm = BoundaryNorm(bounds, cmap.N)



    for label, xaxis, yaxis, xerror, yerror, axis in plot_configs:
        sc = axis.scatter(xaxis, yaxis, c=crossings_per_type, cmap=cmap, norm=norm, s=0.1)


        if label =="xy":
            fig.colorbar(sc, ax=ax, label=r'Number of crossings', boundaries=bounds, spacing='proportional')


        if label=="xy":
            colors = sc.to_rgba(crossings_per_type)


        for x, y, xerr, yerr, color, i in zip(xaxis, yaxis, xerror, yerror, colors, range(len(crossings_per_type))):
            axis.errorbar(x, y, fmt='none', xerr=xerr, yerr=yerr, color=color, marker=None, zorder=crossings_per_type[i])

        
        if label == "Cylindrical":
            plot_magnetospheric_boundaries(axis, frame='MSM', sub_solar_magnetopause=1.4, cylindrical=True, add_legend=True)
            axis.set_xlim(-8,6)
            axis.set_ylim(0,10)
            axis.legend()
        else:
            plot_magnetospheric_boundaries(axis, frame='MSM', sub_solar_magnetopause=1.4, plane=label, add_legend=False)
            axis.set_xlim(-8,8)
            axis.set_ylim(-8,8)


        axis.set_xlabel(f"{xaxis.name} Mercury Radii")
        axis.set_ylabel(f"{yaxis.name} Mercury Radii")


    fig.suptitle(title + " Encounters in MESSENGER Mission")

crossing_per_encounter_config = [
        ("Hollman", hollman_encounters_data),
        ]

for label, data in crossing_per_encounter_config:
    mp_encounters = [i for i in data if "MP" in i['Label'][0]]
    bs_encounters = [i for i in data if "BS" in i['Label'][0]]

    number_of_crossing_encounter = [len(i) for i in data]
    number_of_crossing_encounter_mp = [len(i) for i in mp_encounters]
    number_of_crossing_encounter_bs = [len(i) for i in bs_encounters]

    encounter_bins = np.arange(0, 35, 1)

    configs = [
            ("BS", number_of_crossing_encounter_bs, "yellow"),
            ("MP", number_of_crossing_encounter_mp, "purple"),
            ("total", number_of_crossing_encounter, "red"),
            ]

    encounters_hist = []

    for encounter_type, data_type, color in configs:
        hist = HistogramPanel(data_type, bins=encounter_bins, color=color, minmax=True)

        hist.ax_set_params = {
                "title": f"Number of {encounter_type} crossings per encounter {label}",
                "xlabel": "Number of crossings",
                "ylabel": "Number of encounters",
                "yscale": "log"
                }
        encounters_hist.append(hist)

    encounters_sub_hist = encounters_hist[0] + encounters_hist[1]

    encounters_hist[-1].plot(show=False)
    encounters_sub_hist.plot(show=False)

plt.show()
