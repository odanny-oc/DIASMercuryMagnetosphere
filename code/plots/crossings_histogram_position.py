import numpy as np
import matplotlib.pyplot as plt
from astropy.table import QTable
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

import matplotlib

matplotlib.use('QtAgg')

home_dir = os.getenv('HOME')
data_dir = os.path.join(home_dir, ".ephemeris_data/")

peak_data = QTable.read(data_dir + "peaks_data.csv")

peak_times = Time(peak_data["UTC"]).to_datetime()

print(len(peak_times))

crossing_data = QTable.read(data_dir + "hollman_2025_crossing_list.csv")

crossing_numbers = []

orbit_times =Time(crossing_data["Time"]).to_datetime()

orbit_list = []

delta_t_between_orbits = []


for i in range(len(peak_times) - 1):
    time_start =  peak_times[i]
    time_end = peak_times[i + 1]
    orbit_time = time_end - time_start 
    delta_t_between_orbits.append(orbit_time)
    if orbit_time > timedelta(hours=13) or orbit_time < timedelta(hours=7):
        continue
    else:
        mask = (orbit_times >= peak_times[i]) & (orbit_times <= peak_times[i + 1])
        orbit_list.append(orbit_times[mask].tolist())


for i in range(len(orbit_list)):
    crossing_numbers.append(len(orbit_list[i]))


print(len(orbit_list), len(crossing_data))

average = np.mean(crossing_numbers)

fig, ax = plt.subplots()
ax.tick_params(labelsize=12)
# Plot histogram of number of crosses per orbit

bins = np.arange(0,50,1)
ax.hist(crossing_numbers, bins = bins, edgecolor= 'k')

ax.axvline(average, ls='--', color='r', label=f'Average number of crossings = {average:.2f}')

ax.set_title("Number of magnetic field boundary crossings per orbit, (2011-02-23 - 2015-04-30)", fontsize=24)
ax.set_xlabel("Number of crossings", fontsize=16)
ax.set_ylabel("Number of orbits", fontsize=16)

plot_handles, plot_labels = ax.get_legend_handles_labels()
handles = [Line2D([], [], color = 'none', label = f"Number of orbits = {len(crossing_numbers)}")]
ax.legend(handles= plot_handles + handles, fontsize=14)

print('Plotted Histogram')
plt.show()
