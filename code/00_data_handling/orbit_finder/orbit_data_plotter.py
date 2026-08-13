import os
from pathlib import Path
os.chdir(Path(__file__).resolve().parent)

import matplotlib.pyplot as plt
import numpy as np

import os

from astropy.time import Time
import datetime as dt

from hermpymod.functions.ephemeris_downsampler import parse_crossing_list, parse_periapsis_data
from hermpymod.functions.encounters import parse_encounters_list
from hermpymod.functions.data_per_orbit import orbit_data
from hermpymod.functions.mag_data_plotter import mag_data_plotter
from hermpymod.functions.plot_all import plot_all_ephemeris

home_dir = os.getenv('HOME')
from hermpymod.paths import DATA_DIR
data_dir = DATA_DIR
img_dir = "../../../plots_and_images/"

plt.style.use(img_dir + "presentation.mplstyle")

crossing_data = parse_crossing_list()

orbit_data = orbit_data(crossing_data)

peaks_data = parse_periapsis_data()
peak_times = Time(peaks_data["UTC"]).to_datetime()

orbit_list = [[peak_times[i], peak_times[i+1]] for i in range(len(peak_times) - 1)]

bins = np.arange(0,7, 1)

"""
Get times for each orbit to do the same for Philpott and Sun
"""

hollman_encounters_list = parse_encounters_list()

hollman_time_start = Time(hollman_encounters_list["Time Start"]).to_datetime() 
hollman_time_end=Time(hollman_encounters_list["Time End"]).to_datetime() 


hollman_encounters_per_orbit = []

encounter_lists=[
        (hollman_encounters_per_orbit, hollman_encounters_list, hollman_time_start, hollman_time_end),
        ]


for orbit in orbit_list:
    for encounter_list, data, t_start, t_end in encounter_lists:
        mask_start = (t_start >= orbit[0]) & (t_start <= orbit[-1])
        mask_end = (t_end >= orbit[0]) & (t_end <= orbit[-1])
        mask = mask_start | mask_end

        encounter_list.append(data[mask])


orbit_number = []


for idx, orbit in enumerate(hollman_encounters_per_orbit):
    if len(orbit) < 4:
        orbit_number.append(idx)


print("Number of orbits wiht less than four encounters", len(orbit_number))
slicing = len(orbit_number) // 4

# Ignore first 12 orbits as they are empty due to no MAG data

"""
Plots orbital MAG and ephemeris data
"""

plots = [1103]


for orbit in plots:
    print(f"Plotting orbit number {orbit + 1} ...")
    encounter_times = [orbit_list[orbit][0], orbit_list[orbit][-1]]
    encounter_times = Time(encounter_times).to_datetime()
    tw = dt.timedelta(hours=2)
    time = [encounter_times[0] - tw, encounter_times[-1] + tw]

    fig, ax_mag = mag_data_plotter(time)
    for ax in ax_mag:
        ax.axvspan(orbit_list[orbit][0], orbit_list[orbit][-1], alpha=0.5, color='green', label=f"Orbit number {orbit + 1}")
        ax.legend()

    plot_all_ephemeris(time, crossings=crossing_data, encounters=hollman_encounters_list)
    
    plt.show()
