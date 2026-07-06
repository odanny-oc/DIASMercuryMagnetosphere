import matplotlib.pyplot as plt
import numpy as np

import os

from hermpymod.classes.panels import PlanarplotPanel, HistogramPanel
from hermpymod.functions.ephemeris_downsampler import parse_spice_downsampled, parse_crossing_list
from hermpymod.functions.encounters import encounter_finder

home_dir = os.getenv('HOME')
data_dir = os.path.join(home_dir, '.ephemeris_data/')

plt.style.use("../../plots_and_images/presentation.mplstyle")

crossing_data = parse_crossing_list()

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
