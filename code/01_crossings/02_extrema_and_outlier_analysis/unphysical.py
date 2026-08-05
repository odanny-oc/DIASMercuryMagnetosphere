import os
from pathlib import Path
os.chdir(Path(__file__).resolve().parent)

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.table import QTable, vstack
from sunpy.time import TimeRange
from astropy.time import Time
import pickle

from hermpy.data import parse_messenger_fips, parse_messenger_mag
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
from hermpymod.functions.data_per_orbit import orbit_data
from hermpymod.functions.ephemeris_downsampler import parse_crossing_list, parse_periapsis_data


mpl.use('TkAgg')
home_dir = os.getenv('HOME')
data_dir = os.path.join(home_dir, ".ephemeris_data/")
img_dir = "../../../plots_and_images/"

plt.style.use(img_dir + "presentation.mplstyle")

crossing_data = parse_crossing_list()
crossing_data["UTC"] = Time(crossing_data["UTC"]).to_datetime()

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

unphysical_orbit_crossings = []
unphysical_indices = []

for idx, orbit in enumerate(crossing_orbit_list):
    for label in orbit["Label"]:
        if "UNPHYSICAL" in label:
            unphysical_orbit_crossings.append(orbit)
            unphysical_indices.append(idx)
        else:
            continue


tw = dt.timedelta(hours=3)

print("Number of 'unphysical' crossings ", len(unphysical_orbit_crossings))


"""
Plot all MAG data of orbits with crossings labelled 'unphysical'
"""

for idx, i in enumerate(unphysical_indices):
    times = [crossing_orbit_times[i][0] - tw,crossing_orbit_times[i][-1] + tw]
    print(idx + 1, times)
    fig, ax = mag_data_plotter(times)

    for ax in ax:
        ax.axvspan(orbit_times[i][0], orbit_times[i][-1], alpha=0.3, color= "green", label=f"Orbit {i}\n {idx+1}/{len(unphysical_indices)}")
        ax.legend(fontsize=16)

    plt.show()

