import os
from pathlib import Path
os.chdir(Path(__file__).resolve().parent)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.table import QTable, vstack
from sunpy.time import TimeRange
from astropy.time import Time
import pickle

from hermpy.net import ClientMESSENGER
from hermpy.plotting import MultiPanel, SpectrogramPanel, TimeseriesPanel

import matplotlib.pyplot as plt
from matplotlib import colormaps as cm
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from PIL import Image
import datetime as dt
from datetime import timedelta
import os

from hermpymod.classes.panels import HistogramPanel
from hermpymod.functions.mag_data_plotter import mag_data_plotter
from hermpymod.functions.ephemeris_downsampler import parse_crossing_list, parse_periapsis_data
from hermpymod.functions.data_per_orbit import orbit_data


mpl.use('TkAgg')
home_dir = os.getenv('HOME')
from hermpymod.paths import DATA_DIR
data_dir = DATA_DIR
img_dir = "../../../plots_and_images/"

plt.style.use(img_dir + "presentation.mplstyle")

crossing_data = parse_crossing_list()
crossing_data["UTC"] = Time(crossing_data["UTC"]).to_datetime()

mask = ["UNPHYSICAL" not in label for label in crossing_data["Label"]]
crossing_data = crossing_data[mask]

peaks_data = parse_periapsis_data()
peak_times = Time(peaks_data["UTC"]).to_datetime()

crossing_numbers = []

crossing_times =Time(crossing_data["UTC"]).to_datetime()

crossing_orbit_list = orbit_data(crossing_data=crossing_data)
crossing_orbit_times = [time["UTC"] for time in crossing_orbit_list]

# Start from Hollman list (orbit 12; index 11)
orbit_times = [(peak_times[i], peak_times[i+1]) for i in range(11, len(peaks_data) - 1)]


for i in range(len(crossing_orbit_list)):
    crossing_numbers.append(len(crossing_orbit_list[i]))


crossing_numbers = np.array(crossing_numbers)

types_config = [
    ("BS", ["BS" in crossing_data["Label"][i] for i in range(len(crossing_data))], "yellow"),
    ("MP", ["MP" in crossing_data["Label"][i] for i in range(len(crossing_data))], "purple"),
    ]

bs_mp_orbit_crossings_mask = []


for label, mask, color in types_config:
    time_type = crossing_times[mask]

    num_crossing = np.array([len([time for time in orbit if time in time_type]) for orbit in crossing_orbit_times])

    mask = (num_crossing % 2 == 1)
    indices = [i for i, val in enumerate(mask) if val]
    bs_mp_orbit_crossings_mask.append(indices)


crossing_orbit_times = np.array(crossing_orbit_times, dtype=object)

odd_crossings_mask = (crossing_numbers %2 == 1)
odd_orbit_indices = [i for i, val in enumerate(odd_crossings_mask) if val]
print("Orbit number of orbits with an odd number of crossings ", odd_orbit_indices)

odd_orbit_crossings = [crossing_orbit_list[i] for i in odd_orbit_indices]

bs_odd_orbit_crossings = [[i for i in label["Label"] if i[:2] == "BS"] for label in odd_orbit_crossings]
bs_odd_orbit_crossings = [i for i in bs_odd_orbit_crossings if len(i) %2==1]
bs_odd_orbit_crossings_tot = crossing_orbit_times[bs_mp_orbit_crossings_mask[0]]

mp_odd_orbit_crossings = [[i for i in label["Label"] if i[:2] == "MP"] for label in odd_orbit_crossings]
mp_odd_orbit_crossings = [i for i in mp_odd_orbit_crossings if len(i) %2==1]
mp_odd_orbit_crossings_tot = crossing_orbit_times[bs_mp_orbit_crossings_mask[1]]

print("Number of orbits with odd number of crossings", len(odd_orbit_crossings))

print("Number of BS orbits with odd number of crossings", len(bs_odd_orbit_crossings))

print("Number of MP orbits with odd number of crossings", len(mp_odd_orbit_crossings))

print("Number of BS orbits with odd number of crossings, total", len(bs_odd_orbit_crossings_tot))

print("Number of MP orbits with odd number of crossings, total", len(mp_odd_orbit_crossings_tot))


tw = dt.timedelta(hours=3)


"""
Plots all odd crossings MAG data to verify that they are due to either a data gap or a misplaced/missing crossing.
"""

for idx, i in enumerate(odd_orbit_indices):
    times = [crossing_orbit_times[i][0] - tw,crossing_orbit_times[i][-1] + tw]
    print(idx, times)
    fig, ax = mag_data_plotter(times)

    for ax in ax:
        ax.axvspan(orbit_times[i][0], orbit_times[i][-1], alpha=0.3, color= "green", label=f"Orbit {i}\n {idx+1}/{len(odd_orbit_indices)}")
        ax.legend()

    plt.show()
