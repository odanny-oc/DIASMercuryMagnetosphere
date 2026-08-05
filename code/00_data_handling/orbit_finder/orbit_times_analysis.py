import os
from pathlib import Path
os.chdir(Path(__file__).resolve().parent)

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt

from hermpy.plotting import TimeseriesPanel

from astropy.table import QTable
from sunpy.time import TimeRange
from astropy.time import Time

import datetime as dt
import os

from hermpymod.classes.panels import HistogramPanel, PlanarplotPanel
from hermpymod.functions.ephemeris_downsampler import parse_periapsis_data, parse_crossing_list
from hermpymod.functions.downsampled_positional_data import parse_spice_downsampled

home_dir = os.getenv('HOME')
data_dir = os.path.join(home_dir, ".ephemeris_data/")
img_dir = "../../../plots_and_images/"

mpl.use('QtAgg')

plt.style.use(img_dir + "presentation.mplstyle")

# Load periapsis data
peak_data = parse_periapsis_data()

peak_times = Time(peak_data["UTC"]).to_datetime()

print("Number of peaks", len(peak_times))

# Load crossing data
crossing_data = parse_crossing_list()

crossing_times =Time(crossing_data["UTC"]).to_datetime()

orbit_list = []


"""
Organise crossings by what orbit they are in
"""

for i in range(len(peak_times) - 1):
    time_start =  peak_times[i]
    time_end = peak_times[i + 1]
    orbit_time = time_end - time_start 
    if orbit_time > dt.timedelta(hours=13) or orbit_time < dt.timedelta(hours=7):
        continue
    else:
        mask = (crossing_times >= peak_times[i]) & (crossing_times <= peak_times[i + 1])
        orbit_list.append(crossing_times[mask])


# 10 minute long bins
time_bins = np.arange(6,13, 0.167)

delta_t_numeric_full = peak_data["Orbit Length"]

# Excludes first point which is 0 for array size
delta_t_numeric = np.array(delta_t_numeric_full[1:])

delta_t_mask = (delta_t_numeric <= 7) & (delta_t_numeric >= 13)

print("Length of last orbit in hours (partial)", delta_t_numeric[-1])
print("Number of orbits outside 7-13 hours", len(delta_t_numeric[delta_t_mask]))

transition_orbit_mask = (delta_t_numeric >= 9) & (delta_t_numeric <= 10)

indices = [i for i, val in enumerate(transition_orbit_mask) if val]

print("Orbit number of 9 hour orbits", indices)
print("Number of 9 hour orbits", len(indices))

orbit_list = np.array(orbit_list, dtype=object)
transition_orbits = orbit_list[transition_orbit_mask]

# Start and end times of each 9hr orbit
transition_orbit_times = [(i[0], i[-1]) for i in transition_orbits]

# Times to plot
plotting_times = [transition_orbit_times[0][0]-dt.timedelta(hours=12), transition_orbit_times[-1][-1] + dt.timedelta(hours=12)]

# Parse the orbital data in that time range
orbit_data = parse_spice_downsampled([plotting_times[0], plotting_times[-1]])

time_series_orbit_data = orbit_data["UTC", "|R|"]

time_series_plot = TimeseriesPanel(time_series_orbit_data)

time_series_plot.plot(show=False)


# Plot histogram of orbit times, should be double peaked graph, one for 12hr orbits and one for 8hr orbits
delta_t_between_orbits_hist = HistogramPanel(delta_t_numeric, bins=time_bins, minmax=True)

delta_t_between_orbits_hist.ax_set_params = {
        "title": f"Length of all Orbits ({crossing_times[0]} - {crossing_times[-1]})",
        "xlabel": "Orbit length (hours)",
        "ylabel": "Number of orbits",
        "yscale": "log",
        }


delta_t_between_orbits_hist.plot(show=False)
plt.savefig(img_dir + "hollman_orbit_times.svg")
plt.show()
