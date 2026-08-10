import os
from pathlib import Path
os.chdir(Path(__file__).resolve().parent)

import matplotlib.pyplot as plt
import numpy as np

import os

from astropy.table import QTable
from astropy.time import Time

from hermpymod.classes.panels import PlanarplotPanel, HistogramPanel
from hermpymod.functions.ephemeris_downsampler import parse_crossing_list
from hermpymod.functions.encounters import encounter_finder, parse_encounters_list
from hermpymod.functions.data_per_orbit import orbit_data
from hermpymod.functions.ephemeris_downsampler import parse_periapsis_data

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

"""
Get times for each orbit to do the same for Philpott and Sun
"""

hollman_encounters_data = encounter_finder(crossing_data)
print(len(hollman_encounters_data))
hollman_encounters_list = parse_encounters_list()
print(len(hollman_encounters_list))

hollman_time_start = Time(hollman_encounters_list["Time Start"]).to_datetime() 
hollman_time_end=Time(hollman_encounters_list["Time End"]).to_datetime() 

sun_encounters_data = QTable.read(data_dir + "sun_2023_crossing.csv")
sun_time_start = Time(sun_encounters_data["Time Start"]).to_datetime() 
sun_time_end=Time(sun_encounters_data["Time End"]).to_datetime() 

philpott_encounters_data = QTable.read(data_dir + "philpott_encounter_list_2020.csv")
phil_time_start = Time(philpott_encounters_data["Time Start"]).to_datetime()
phil_time_end = Time(philpott_encounters_data["Time End"]).to_datetime()


hollman_encounters_per_orbit = []
sun_encounters_per_orbit = []
phil_encounters_per_orbit = []

encounter_lists=[
        (hollman_encounters_per_orbit, hollman_encounters_list, hollman_time_start, hollman_time_end),
        (sun_encounters_per_orbit, sun_encounters_data, sun_time_start, sun_time_end),
        (phil_encounters_per_orbit, philpott_encounters_data, phil_time_start, phil_time_end),
        ]


for orbit in orbit_list:
    for encounter_list, data, t_start, t_end in encounter_lists:
        mask_start = (t_start >= orbit[0]) & (t_start <= orbit[-1])
        mask_end = (t_end >= orbit[0]) & (t_end <= orbit[-1])
        mask = mask_start | mask_end

        encounter_list.append(data[mask])

bins = np.arange(0,7, 1)


encounters_per_orbit = [
        ("Hollman", hollman_encounters_per_orbit, hollman_encounters_list),
        ("Sun et al.", sun_encounters_per_orbit, sun_encounters_data),
        ("Philpott et al.", phil_encounters_per_orbit, philpott_encounters_data)
        ]


for label, data, encounters in encounters_per_orbit:

    """
    Remove beginning orbits with no data
    """
    i = 0
    while i < len(data) and len(data[i]) == 0:
        i += 1
    data = data[i:]

    num_per_orbit = [len(i) for i in data]

    overcount = sum(num_per_orbit) - len(encounters)
    print(label + " Overcounted encounters", sum(num_per_orbit), " - ", len(encounters), " = ", overcount)
    encounters_per_orbit_hist = HistogramPanel(num_per_orbit, bins)
    encounters_per_orbit_hist.ax_set_params = {
            "title": f"Number of encounters per orbit {label}",
            "xlabel": "Number of encounters",
            "ylabel": "Number of orbits",
            "yscale": "log"
            }
    encounters_per_orbit_hist.plot(show=False)
    plt.savefig(img_dir + f"label_encounters_per_orbit.svg")

plt.show()
