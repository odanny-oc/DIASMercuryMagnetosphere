import os
from pathlib import Path
os.chdir(Path(__file__).resolve().parent)

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.gridspec import GridSpec

import matplotlib as mpl
import matplotlib.pyplot as plt

from hermpy.plotting import MultiPanel, SpectrogramPanel, TimeseriesPanel

from astropy.table import QTable
from sunpy.time import TimeRange
from astropy.time import Time

import datetime as dt
import os

from hermpymod.classes.panels import HistogramPanel, PlanarplotPanel
from hermpymod.functions.mag_data_plotter import mag_data_plotter
from hermpymod.functions.data_per_orbit import orbit_data
from hermpymod.functions.plot_all import plot_all_ephemeris
from hermpymod.functions.ephemeris_downsampler import parse_crossing_list, parse_periapsis_data
from hermpymod.functions.encounters import parse_encounters_list

home_dir = os.getenv('HOME')
from hermpymod.paths import DATA_DIR
data_dir = DATA_DIR
img_dir = "../../../plots_and_images/"

mpl.use('QtAgg')

plt.style.use(img_dir + "presentation.mplstyle")

peak_data = parse_periapsis_data()

peak_times = Time(peak_data["UTC"]).to_datetime()

# Start from Hollman list (orbit 12; index 11)
orbit_times = [(peak_times[i], peak_times[i+1]) for i in range(11, len(peak_data) - 1)]

print("Number of peaks found", len(peak_times))

crossing_data = parse_crossing_list()
crossing_data["UTC"] = Time(crossing_data["UTC"]).to_datetime()

encounters_data = parse_encounters_list()
encounters_data["Time Start"] = Time(encounters_data["Time Start"]).to_datetime()
encounters_data["Time End"] = Time(encounters_data["Time End"]).to_datetime()

crossing_times =Time(crossing_data["UTC"]).to_datetime()

crossing_orbit_list = orbit_data(crossing_data=crossing_data)
crossing_orbit_times = [orbit["UTC"] for orbit in crossing_orbit_list]

types_config = [
    ("BS", ["BS" in crossing_data["Label"][i] for i in range(len(crossing_data))], "yellow"),
    # ("MP", ["MP" in crossing_data["Label"][i] for i in range(len(crossing_data))], "purple"),
        ]

encounters_per_orbit = []


for orbit in orbit_times:
    mask_start = (encounters_data["Time Start"] >= orbit[0]) & (encounters_data["Time Start"] <= orbit[-1])
    mask_end = (encounters_data["Time End"] >= orbit[0]) & (encounters_data["Time End"] <= orbit[-1])
    mask = mask_start | mask_end

    encounters_per_orbit.append(encounters_data[mask])


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

        fig, ax = plot_all_ephemeris(times[i], crossings=crossing_data, encounters=encounters_data)

        for num in range(3):
            ax[num].set_xlim(-5,5)
            ax[num].set_ylim(-6,4)

        fig.suptitle(f"{times[i][0].isoformat()}-{times[i][-1].isoformat()} Number of crossings {max_list[i]}", fontsize=20)

        plt.savefig(img_dir + f"orbit_with_{max_list[i]}_crossings.svg")

        tw = dt.timedelta(hours=3)
        mag_times = (orbit_times[max_list_indices[i]][0] - tw, orbit_times[max_list_indices[i]][-1] + tw)
        mag_times_zoom = (encounters_per_orbit[max_list_indices[i]]["Time Start"][1], encounters_per_orbit[max_list_indices[i]]["Time End"][-2])

        fig, ax = mag_data_plotter(mag_times, zoom=mag_times_zoom)


        for ax in ax:
            ax.axvspan(orbit_times[max_list_indices[i]][0], orbit_times[max_list_indices[i]][-1], alpha=0.3, color= "green", label=f"Orbit {max_list_indices[i]}\n Number of crossings {max_list[i]}")
            ax.legend(fontsize=16, loc="center left")


plt.show()
