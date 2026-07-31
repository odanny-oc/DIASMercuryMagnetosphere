import matplotlib.pyplot as plt
import numpy as np

import os

from astropy.table import QTable
from astropy.time import Time
import datetime as dt

from hermpymod.classes.panels import PlanarplotPanel, HistogramPanel
from hermpymod.functions.ephemeris_downsampler import parse_crossing_list, parse_periapsis_data
from hermpymod.functions.encounters import encounter_finder, parse_encounters_list
from hermpymod.functions.data_per_orbit import orbit_data
from hermpymod.functions.mag_data_plotter import mag_data_plotter
from hermpymod.functions.plot_all import plot_all_ephemeris

home_dir = os.getenv('HOME')
data_dir = os.path.join(home_dir, '.ephemeris_data/')
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

phil_encounters_list = QTable.read(data_dir + "philpott_encounter_list_2020.csv")

phil_time_start = Time(phil_encounters_list["Time Start"]).to_datetime() 
phil_time_end=Time(phil_encounters_list["Time End"]).to_datetime() 


phil_encounters_per_orbit = []

encounter_lists=[
        (phil_encounters_per_orbit, phil_encounters_list, phil_time_start, phil_time_end),
        ]


for orbit in orbit_list:
    for encounter_list, data, t_start, t_end in encounter_lists:
        mask_start = (t_start >= orbit[0]) & (t_start <= orbit[-1])
        mask_end = (t_end >= orbit[0]) & (t_end <= orbit[-1])
        mask = mask_start | mask_end

        encounter_list.append(data[mask])

overcounted = []
orbit_number =[]

for idx, orbit in enumerate(phil_encounters_per_orbit):
    try:
        if orbit[-1] == phil_encounters_per_orbit[idx + 1][0]:
            overcounted.append(orbit[-1])
            orbit_number.append(idx)
    # Skip orbits with no encounters and final orbit
    except IndexError:
        continue

print(len(overcounted))


for encounter, orbit in zip(overcounted, orbit_number):
    encounter_times = [encounter["Time Start"], encounter["Time End"]]
    encounter_times = Time(encounter_times).to_datetime()
    tw = dt.timedelta(hours=4)
    time = [encounter_times[0] - tw, encounter_times[-1] + tw]

    fig, ax_mag = mag_data_plotter(time)
    for ax in ax_mag:
        ax.axvspan(encounter_times[0], encounter_times[-1], alpha=0.5, color='orange', label="Encounter")
        ax.axvspan(orbit_list[orbit][0], orbit_list[orbit][-1], alpha=0.5, color='green', label=f"Orbit number {orbit + 1}")

    plot_all_ephemeris(time)
    
plt.show()
