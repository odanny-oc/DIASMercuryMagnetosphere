import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt

from hermpy.plotting import MultiPanel, SpectrogramPanel, TimeseriesPanel

from astropy.table import QTable
from sunpy.time import TimeRange
from astropy.time import Time

import datetime as dt
import os

from hermpymod.classes.panels import HistogramPanel, PlanarplotPanel
from hermpymod.functions.ephemeris_downsampler import parse_spice_downsampled
from hermpymod.functions.mag_data_plotter import mag_data_plotter
from hermpymod.functions.data_per_orbit import orbit_data

home_dir = os.getenv('HOME')
data_dir = os.path.join(home_dir, ".ephemeris_data/")
img_dir = "../../plots_and_images/"

mpl.use('QtAgg')

plt.style.use(img_dir + "presentation.mplstyle")

peak_data = QTable.read(data_dir + "peaks_data.csv")

peak_times = Time(peak_data["UTC"]).to_datetime()

# Start from Hollman list (orbit 12; index 11)
orbit_times = [(peak_times[i], peak_times[i+1]) for i in range(11, len(peak_data) - 1)]

print("Number of peaks found", len(peak_times))

crossing_data = QTable.read(data_dir + "hollman_2025_crossing_list.ecsv")

crossing_times =Time(crossing_data["UTC"]).to_datetime()

crossing_orbit_list = orbit_data(crossing_data=crossing_data)
crossing_orbit_times = [orbit["UTC"].to_datetime() for orbit in crossing_orbit_list]

types_config = [
    ("BS", ["BS" in crossing_data["Label"][i] for i in range(len(crossing_data))], "yellow"),
    # ("MP", ["MP" in crossing_data["Label"][i] for i in range(len(crossing_data))], "purple"),
        ]

total_crossing_numbers = [len(i) for i in crossing_orbit_list]

histograms = []

for label, mask, color in types_config:
    time_type = crossing_times[mask]

    num_crossing_from_orbit = np.array([len([time for time in orbit if time in time_type]) for orbit in crossing_orbit_times])
    max_index = np.where(num_crossing_from_orbit==np.max(num_crossing_from_orbit))[0][0]
    second_largest = sorted(num_crossing_from_orbit)[-2]
    second_max_index = np.where(num_crossing_from_orbit==second_largest)[0][0]

    max_list = [num_crossing_from_orbit[max_index], second_largest]
    max_list_indices = [max_index, second_max_index]

    times = [(crossing_orbit_times[max_index][0], crossing_orbit_times[max_index][-1]), (crossing_orbit_times[second_max_index][0], crossing_orbit_times[second_max_index][-1])]
    for i in range(len(times)):
        planar_plot_xy = PlanarplotPanel(times[i], "X-Y", crossings=True)
        planar_plot_xz = PlanarplotPanel(times[i], "X-Z", crossings=True)
        plots = planar_plot_xy + planar_plot_xz
        plots.ax_set_param = {
                "title": f"{times[i]} Number of crossings {max_list[i]}" 
                }
        plots.plot(show=False)

        tw = dt.timedelta(hours=3)
        mag_times = (orbit_times[max_list_indices[i]][0] - tw, orbit_times[max_list_indices[i]][-1] + tw)
        fig, ax = mag_data_plotter(mag_times)

        for ax in ax:
            ax.axvspan(orbit_times[max_list_indices[i]][0], orbit_times[max_list_indices[i]][-1], alpha=0.3, color= "green", label=f"Orbit {max_list_indices[i]}\n Number of crossings {max_list[i]}")
            ax.legend(fontsize=16)


plt.show()
