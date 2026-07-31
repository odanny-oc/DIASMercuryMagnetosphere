import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.table import QTable
from sunpy.time import TimeRange
from astropy.time import Time

from hermpy.plotting import MultiPanel, SpectrogramPanel, TimeseriesPanel

import matplotlib.pyplot as plt
from datetime import timedelta
import os

from hermpymod.classes.panels import HistogramPanel
from hermpymod.functions.data_per_orbit import orbit_data
from hermpymod.functions.ephemeris_downsampler import parse_periapsis_data

mpl.use('QtAgg')

home_dir = os.getenv('HOME')
data_dir = os.path.join(home_dir, ".ephemeris_data/")
img_dir = "../../../plots_and_images/"

plt.style.use(img_dir + "presentation.mplstyle")

peaks_data = parse_periapsis_data()

peak_times = Time(peaks_data["UTC"]).to_datetime()

print("Number of periapsides", len(peak_times))

crossing_data = QTable.read(data_dir + "hollman_2025_crossing_list.ecsv")

crossing_numbers = []

orbit_times =Time(crossing_data["UTC"]).to_datetime()

orbit_list = []

delta_t_between_orbits = []


orbit_list = orbit_data()

for i in range(len(orbit_list)):
    crossing_numbers.append(len(orbit_list[i]))


print(len(orbit_list), len(crossing_data))

average = np.mean(crossing_numbers)

# Plot histogram of number of crosses per orbit

bins = np.arange(0,50,1)

crossings_per_orbit_hist = HistogramPanel(crossing_numbers, bins=bins)

crossings_per_orbit_hist.ax_set_params = {
        "title": f"Number of magnetic field boundary crossings per orbit, ({Time(peak_times[0]).iso[:11]} - {Time(peak_times[-1]).iso[:11]})",
        "xlabel": "Number of crossings",
        "ylabel":"Number of orbits",
        "yscale":"log",
}

print('Plotted Histogram')


crossings_per_orbit_hist.plot(show=False)
# plt.savefig(img_dir + "crossings_per_orbit_hollman.svg")
plt.savefig(img_dir + "crossings_per_orbit_hollman_log.svg")
plt.show()
