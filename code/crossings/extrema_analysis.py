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
        planar_plot_xy = PlanarplotPanel(times[i], "X-Y", crossings=True, encounters=True)
        planar_plot_xz = PlanarplotPanel(times[i], "X-Z", crossings=True, encounters=True)
        planar_plot_yz = PlanarplotPanel(times[i], "Y-Z", crossings=True, encounters=True)
        planar_plot_cylind = PlanarplotPanel(times[i], "Y-Z", cylindrical=True, crossings=True, encounters=False)
        plots_list = [planar_plot_xy , planar_plot_xz , planar_plot_yz , planar_plot_cylind]
        plots = planar_plot_xy + planar_plot_xz + planar_plot_yz + planar_plot_cylind
        fig = plt.figure()
        gs = GridSpec(3, 2, figure=fig)

        ax = [fig.add_subplot(gs[0, 0]),
         fig.add_subplot(gs[0, 1]),
         fig.add_subplot(gs[1, :]), 
         fig.add_subplot(gs[2, :]),]
        # fig, ax = plots.plot(show=False,sharex=False)
        for plot, axis in zip(plots_list, ax):
            plot._plot_on(axis)

        handles, _ = ax[0].get_legend_handles_labels()

        for num in range(3):
            ax[num].set_xlim(-5,5)
            ax[num].set_ylim(-6,4)

        for ax in ax:
            ax.xaxis.label.set_size(14)
            ax.yaxis.label.set_size(14)
            ax.tick_params(axis='both', labelsize=12)
            ax.legend().remove()



        fig.legend(handles=handles, bbox_to_anchor=(1.0,0.9), fontsize=11)
        fig.suptitle(f"{times[i][0].isoformat()}-{times[i][-1].isoformat()} Number of crossings {max_list[i]}", fontsize=20)

        plt.savefig(img_dir + f"orbit_with_{max_list[i]}_crossings.svg")

        tw = dt.timedelta(hours=3)
        mag_times = (orbit_times[max_list_indices[i]][0] - tw, orbit_times[max_list_indices[i]][-1] + tw)
        fig, ax = mag_data_plotter(mag_times)


        for ax in ax:
            ax.axvspan(orbit_times[max_list_indices[i]][0], orbit_times[max_list_indices[i]][-1], alpha=0.3, color= "green", label=f"Orbit {max_list_indices[i]}\n Number of crossings {max_list[i]}")
            ax.legend(fontsize=16)


plt.show()
