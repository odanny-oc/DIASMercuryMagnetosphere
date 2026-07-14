import matplotlib.pyplot as plt
import numpy as np

import os

from hermpymod.classes.panels import PlanarplotPanel, HistogramPanel
from hermpymod.functions.ephemeris_downsampler import parse_crossing_list
from hermpymod.functions.encounters import encounter_finder
from hermpymod.functions.data_per_orbit import orbit_data

home_dir = os.getenv('HOME')
data_dir = os.path.join(home_dir, '.ephemeris_data/')
img_dir = "../../plots_and_images/"

plt.style.use("../../plots_and_images/presentation.mplstyle")

crossing_data = parse_crossing_list()

orbit_data = orbit_data(crossing_data)

encounters_per_orbit = [encounter_finder(orbit) for orbit in orbit_data]

encounters_per_orbit = np.array(encounters_per_orbit, dtype=object)

num_encounters_per_orbit = [len(encounter_orbit) for encounter_orbit in encounters_per_orbit]

print(sum(num_encounters_per_orbit))

bins = np.arange(0,7, 1)
encounters_per_orbit_hist = HistogramPanel(num_encounters_per_orbit, bins)
encounters_per_orbit_hist.ax_set_params = {
        "title": f"Number of encounters per orbit",
        "xlabel": "Number of encounters",
        "ylabel": "Number of orbits",
        "yscale": "log"
        }
encounters_per_orbit_hist.plot(show=False)

plt.savefig(img_dir + "encounters_per_orbit.svg")

encounters_data = encounter_finder(crossing_data)

mp_encounters = [i for i in encounters_data if "MP" in i['Label'][0]]
bs_encounters = [i for i in encounters_data if "BS" in i['Label'][0]]

number_of_crossing_encounter = [len(i) for i in encounters_data]
number_of_crossing_encounter_mp = [len(i) for i in mp_encounters]
number_of_crossing_encounter_bs = [len(i) for i in bs_encounters]

encounter_bins = np.arange(0, 17, 1)

configs = [
        ("BS", number_of_crossing_encounter_bs, "yellow"),
        ("MP", number_of_crossing_encounter_mp, "purple"),
        ("total", number_of_crossing_encounter, "red"),
        ]

encounters_hist = []

for label, data, color in configs:
    hist = HistogramPanel(data, bins=encounter_bins, color=color)

    hist.ax_set_params = {
            "title": f"Number of {label} crossings per encounter (defined by crossings of the same type and direction)",
            "xlabel": "Number of crossings",
            "ylabel": "Number of encounters",
            "yscale": "log"
            }
    encounters_hist.append(hist)

encounters_sub_hist = encounters_hist[0] + encounters_hist[1]

encounters_hist[-1].plot(show=False)
encounters_sub_hist.plot(show=False)

plt.show()
