import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd

import matplotlib as mpl
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from matplotlib import colormaps as cm

from hermpy.plotting import MultiPanel, SpectrogramPanel, TimeseriesPanel

from astropy.time import TimeDelta
from astropy.table import QTable, vstack
from sunpy.time import TimeRange
from astropy.time import Time

import datetime as dt
import os

from hermpymod.classes.panels import HistogramPanel, PlanarplotPanel
from hermpymod.functions.ephemeris_downsampler import parse_spice_downsampled
from hermpymod.functions.data_per_orbit import orbit_data


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

crossings_per_orbit_list = orbit_data()

crossings_per_orbit_times = [orbit["UTC"].to_datetime() for orbit in crossings_per_orbit_list]

delta_t_between_orbits = pd.to_timedelta(peak_data["delta t"])

types_config = [
    ("BS", ["BS" in crossing_data["Label"][i] for i in range(len(crossing_data))], "yellow"),
    ("MP", ["MP" in crossing_data["Label"][i] for i in range(len(crossing_data))], "purple"),
    ("IN", ["IN" in crossing_data["Label"][i] for i in range(len(crossing_data))], "blue"),
    ("OUT", ["OUT" in crossing_data["Label"][i] for i in range(len(crossing_data))], "orange"),
        ]

total_crossing_numbers = [len(i) for i in crossings_per_orbit_list]

histograms = []

for label, mask, color in types_config:
    time_type = crossing_times[mask]

    num_crossing_from_orbit = [len([time for time in orbit if time in time_type]) for orbit in crossings_per_orbit_times]
    hist = HistogramPanel(num_crossing_from_orbit, bins='auto', color=color)

    hist.ax_set_params = {
            "title" : f"Number of {label} crossings per orbit, ({crossing_times[0].isoformat()[:10]} - {crossing_times[-1].isoformat()[:10]})",
            "ylabel": f"Number of orbits",
            "xlabel": f"Number of {label} crossings",
            "yscale": "log",
            }
    histograms.append(hist)

bs_mp_hist = histograms[0] + histograms[1]
in_out_hist = histograms[2] + histograms[3]

bs_mp_hist.plot(show=False, figsize = (18,16))
plt.savefig(img_dir + "hollman_bs_mp_per_orbit.svg")

in_out_hist.plot(show=False)

# Check if all crossings are counted
print("Sum of all crossing counted in the orbit list", sum(total_crossing_numbers))
print("Length of crossing list", len(crossing_data))

plt.show()
