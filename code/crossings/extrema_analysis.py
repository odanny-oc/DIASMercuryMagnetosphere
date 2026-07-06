import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import matplotlib as mpl
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from matplotlib import colormaps as cm
from mpl_toolkits.mplot3d.art3d import Line3DCollection

from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from hermpy.plotting import MultiPanel, SpectrogramPanel, TimeseriesPanel

from astropy.time import TimeDelta
from astropy.table import QTable
from sunpy.time import TimeRange
from astropy.time import Time

import datetime as dt
import sys
import os
import pickle

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from hermpymod.classes.panels import HistogramPanel, PlanarplotPanel
from hermpymod.functions.ephemeris_downsampler import parse_spice_downsampled

home_dir = os.getenv('HOME')
data_dir = os.path.join(home_dir, ".ephemeris_data/")
img_dir = "../../plots_and_images/"

mpl.use('QtAgg')

plt.style.use(img_dir + "presentation.mplstyle")

peak_data = QTable.read(data_dir + "peaks_data.csv")

peak_times = Time(peak_data["UTC"]).to_datetime()

print("Number of peaks found", len(peak_times))

crossing_data = QTable.read(data_dir + "hollman_2025_crossing_list.ecsv")

crossing_times =Time(crossing_data["UTC"]).to_datetime()

orbit_list = []


for i in range(len(peak_times) - 1):
    time_start =  peak_times[i]
    time_end = peak_times[i + 1]
    orbit_time = time_end - time_start 
    if orbit_time > dt.timedelta(hours=13) or orbit_time < dt.timedelta(hours=7):
        continue
    else:
        mask = (crossing_times >= peak_times[i]) & (crossing_times <= peak_times[i + 1])
        orbit_list.append(crossing_times[mask])


delta_t_between_orbits = pd.to_timedelta(peak_data["delta t"])

types_config = [
    ("BS", ["BS" in crossing_data["Label"][i] for i in range(len(crossing_data))], "yellow"),
    # ("MP", ["MP" in crossing_data["Label"][i] for i in range(len(crossing_data))], "purple"),
        ]

total_crossing_numbers = [len(i) for i in orbit_list]

histograms = []

for label, mask, color in types_config:
    time_type = crossing_times[mask]

    num_crossing_from_orbit = np.array([len([time for time in orbit if time in time_type]) for orbit in orbit_list])
    max_index = np.where(num_crossing_from_orbit==np.max(num_crossing_from_orbit))[0][0]
    second_largest = sorted(num_crossing_from_orbit)[-2]
    second_max_index = np.where(num_crossing_from_orbit==second_largest)[0][0]

    max_list = [num_crossing_from_orbit[max_index], second_largest]

    times = [(orbit_list[max_index][0], orbit_list[max_index][-1]), (orbit_list[second_max_index][0], orbit_list[second_max_index][-1])]
    for i in range(len(times)):
        planar_plot_xy = PlanarplotPanel(times[i], "X-Y", crossings=True)
        planar_plot_xz = PlanarplotPanel(times[i], "X-Z", crossings=True)
        plots = planar_plot_xy + planar_plot_xz
        plots.ax_set_param = {
                "title": f"{times[i]} Number of crossings {max_list[i]}" 
                }
        plots.plot(show=False)

plt.show()
