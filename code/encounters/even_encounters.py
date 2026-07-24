import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec


from astropy.table import QTable, vstack, Column
import astropy.units as u
from astropy.time import Time
from astropy.table import QTable, hstack, vstack
from sunpy.time import TimeRange

from hermpy.utils import Constants as c

import datetime as dt
import os

from hermpymod.classes.panels import PlanarplotPanel, HistogramPanel
from hermpymod.functions.ephemeris_downsampler import parse_crossing_list
from hermpymod.functions.encounters import parse_encounters_list, encounter_finder
from hermpymod.functions.crossings_plotter import plot_crossings
from hermpymod.functions.mag_data_plotter import mag_data_plotter
from hermpymod.functions.plot_all import plot_all_ephemeris


home_dir = os.getenv('HOME')
data_dir = os.path.join(home_dir, '.ephemeris_data/')

plt.style.use("../../plots_and_images/presentation.mplstyle")

crossing_data = parse_crossing_list()

crossing_label = crossing_data["Label"]
crossing_direction = crossing_data["Trajectory Direction"]
crossing_time = crossing_data["UTC"].to_datetime()

encounters_data = encounter_finder(crossing_data)
encounters_list = parse_encounters_list()
encounter_start_times = Time(encounters_list["Time Start"]).to_datetime()

# mag_encounters = [i for i in encounters_data if len(i) >= 2]
mag_encounters = [i for i in encounters_data if len(i) %2 == 0]

print(len(mag_encounters))

slicing= len(mag_encounters)//3

mag_encounters = mag_encounters[::slicing]
# mag_encounters.append(encounters_data[i0])

"""
MAG plots for encounters
"""
from astropy.table import Table

for i in mag_encounters:
    Table.pprint(i, max_lines=-1, max_width=-1)
    encounter_time = i["UTC"].to_datetime()
    timewindow = dt.timedelta(minutes=300)
    t_start, t_end = encounter_time[0] - timewindow, encounter_time[-1] + timewindow
    
    i0 = np.searchsorted(encounter_start_times, encounter_time[0], side="left")

    # Plot MAG data in time interval and all crossings
    fig_mag, ax_mag = mag_data_plotter([t_start, t_end])
    print("Mag data plotted")

    # Plot encounter on graph
    ax_mag[0].axvspan(encounter_time[0], encounter_time[-1], alpha=0.5, color='orange', label=f'Encounter number {i0}')
    ax_mag[1].axvspan(encounter_time[0], encounter_time[-1], alpha=0.5, color='orange', label=f'Encounter number {i0}')

    plot_all_ephemeris([t_start, t_end], crossings=crossing_data, encounters=encounters_list)

plt.show()
