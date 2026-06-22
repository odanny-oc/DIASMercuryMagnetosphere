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
import matplotlib
from matplotlib.lines import Line2D
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from hermpymod.classes.panels import HistogramPanel

home_dir = os.getenv('HOME')
data_dir = os.path.join(home_dir, ".ephemeris_data/")

matplotlib.use('QtAgg')

# print(plt.rcParams.keys())

plt.rcParams.update({
    'axes.labelsize' : 15,
    'xtick.labelsize' : 12,
    'ytick.labelsize' : 12,
    'axes.titlesize' : 20,
    'legend.fontsize': 16,
     })

peak_data = QTable.read(data_dir + "peaks_data.csv")

peak_times = Time(peak_data["UTC"]).to_datetime()

print(len(peak_times))

crossing_data = QTable.read(data_dir + "hollman_2025_crossing_list.csv")

bs_mask = ["BS" in crossing_data["Label"][i] for i in range(len(crossing_data))] 
mp_mask = ["MP" in crossing_data["Label"][i] for i in range(len(crossing_data))] 

in_mask = ["IN" in crossing_data["Label"][i] for i in range(len(crossing_data))] 
out_mask = ["OUT" in crossing_data["Label"][i] for i in range(len(crossing_data))] 

crossing_times =Time(crossing_data["Time"]).to_datetime()

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
        mask = (crossing_times >= peak_times[i]) & (crossing_times <= peak_times[i + 1])
        # i0 = np.searchsorted(crossing_times, peak_times[i], side='left')
        # i1 = np.searchsorted(crossing_times, peak_times[i+1], side='right')
        # orbit_list.append(crossing_times[i0:i1])
        orbit_list.append(crossing_times[mask])


time_bins = np.arange(6,13, 0.167)

delta_t_numeric = [dt.total_seconds()/ 3600 for dt in delta_t_between_orbits]

delta_t_between_orbits_hist = HistogramPanel(delta_t_numeric, bins=time_bins)

counts, _ = np.histogram(delta_t_numeric, bins=time_bins)
print(np.min(counts[counts>0]), np.max(counts))

delta_t_between_orbits_hist.ax_set_params = {
        "title": "Orbit length for each orbit",
        "xlabel": "Orbit length (hours)",
        "ylabel": "Number of orbits",
        "yscale": "log",
        }

delta_t_between_orbits_hist.plot(show=False)

total_crossing_numbers = [len(i) for i in orbit_list]

mp_times = crossing_times[mp_mask]

mp_crossings = [len([time for time in orbit if time in mp_times]) for orbit in orbit_list]

mp_hist = HistogramPanel(mp_crossings, bins='auto', color='purple')

bs_times = crossing_times[bs_mask]

bs_crossings = [len([time for time in orbit if time in bs_times]) for orbit in orbit_list]

bs_hist = HistogramPanel(bs_crossings, bins='auto', color='yellow')

mp_hist.ax_set_params = {
        "title" : "Number of magnetopause boundary crossings per orbit, (2011-02-23 - 2015-04-30)",
        "xlabel": "Number of orbits",
        "ylabel": "Number of magnetopause crossings",
        "yscale": "log",
        }

bs_hist.ax_set_params = {
        "title" : "Number of bow shock boundary crossings per orbit, (2011-02-23 - 2015-04-30)",
        "xlabel": "Number of orbits",
        "ylabel": "Number of bow shock crossings",
        "yscale": "log",
        }

bs_mp_hist = bs_hist + mp_hist

in_times = crossing_times[in_mask]

in_crossings = [len([time for time in orbit if time in in_times]) for orbit in orbit_list]

in_hist = HistogramPanel(in_crossings, bins='auto', color='blue')

out_times = crossing_times[out_mask]

out_crossings = [len([time for time in orbit if time in out_times]) for orbit in orbit_list]

out_hist = HistogramPanel(out_crossings, bins='auto', color='orange')

in_hist.ax_set_params = {
        "title" : "Number of IN boundary crossings per orbit, (2011-02-23 - 2015-04-30)",
        "xlabel": "Number of orbits",
        "ylabel": "Number of IN crossings",
        "yscale": "log",
        }

out_hist.ax_set_params = {
        "title" : "Number of OUT boundary crossings per orbit, (2011-02-23 - 2015-04-30)",
        "xlabel": "Number of orbits",
        "ylabel": "Number of OUT crossings",
        "yscale": "log",
        }

in_out_hist = in_hist + out_hist

fig, ax = bs_mp_hist.plot(show=False)

fig, ax = in_out_hist.plot(show=False)

print(sum(total_crossing_numbers), len(crossing_data))

# Plot histogram of number of crosses per orbit

bins = np.arange(0,50,1)
total_hist = HistogramPanel(total_crossing_numbers, bins=bins)

total_hist.ax_set_params = {
        "title" : "Number of magnetic field boundary crossings per orbit, (2011-02-23 - 2015-04-01)",
        "xlabel": "Number of orbits",
        "ylabel": "Number of crossings",
        "yscale": "log",
        }

# ax.set_title("Number of magnetic field boundary crossings per orbit, (2011-02-23 - 2012-03-31)", fontsize=24)
# ax.set_xlabel("Number of crossings", fontsize=16)
# ax.set_ylabel("Number of orbits", fontsize=16)
fig, ax = total_hist.plot(show=False)

plt.show()
