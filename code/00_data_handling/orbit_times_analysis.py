import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt

from hermpy.plotting import TimeseriesPanel

from astropy.table import QTable
from sunpy.time import TimeRange
from astropy.time import Time

import datetime as dt
import os

from hermpymod.classes.panels import HistogramPanel, PlanarplotPanel
from hermpymod.functions.ephemeris_downsampler import parse_periapsis_data
from hermpymod.functions.downsampled_positional_data import parse_spice_downsampled

home_dir = os.getenv('HOME')
data_dir = os.path.join(home_dir, ".ephemeris_data/")
img_dir = "../../plots_and_images/"

mpl.use('QtAgg')

plt.style.use(img_dir + "presentation.mplstyle")

peak_data = parse_periapsis_data()

peak_times = Time(peak_data["UTC"]).to_datetime()

print("Number of peaks", len(peak_times))

crossing_data = QTable.read(data_dir + "hollman_2025_crossing_list.ecsv")

crossing_times =Time(crossing_data["UTC"]).to_datetime()

orbit_list = []


for i in range(len(peak_times) - 1):
    time_start =  peak_times[i]
    time_end = peak_times[i + 1]
    orbit_time = time_end - time_start 
    if orbit_time > dt.timedelta(hours=13) or orbit_time < dt.timedelta(hours=7):
        continue
    else:
        mask = (crossing_times >= peak_times[i]) & (crossing_times <= peak_times[i + 1])
        orbit_list.append(crossing_times[mask])


delta_t_between_orbits = pd.to_timedelta(peak_data["Orbit Length"])


time_bins = np.arange(6,13, 0.167)

delta_t_numeric_full = [dt.total_seconds()/ 3600 for dt in delta_t_between_orbits]

# Excludes first point which is 0 for array size
delta_t_numeric = np.array(delta_t_numeric_full[1:])

delta_t_mask = (delta_t_numeric <= 7) & (delta_t_numeric >= 13)

print("Length of last orbit in hours (partial)", delta_t_numeric[-1])
print("Number of orbits outside 7-13 hours", len(delta_t_numeric[delta_t_mask]))

transition_orbit_mask = (delta_t_numeric >= 9) & (delta_t_numeric <= 10)

indices = [i for i, val in enumerate(transition_orbit_mask) if val]

print("Orbit number of 9 hour orbits", indices)
print("Number of 9 hour orbits", len(indices))

orbit_list = np.array(orbit_list, dtype=object)
transition_orbits = orbit_list[transition_orbit_mask]

transition_orbit_times = [(i[0], i[-1]) for i in transition_orbits]

plotting_times = [
        (transition_orbit_times[0][0]-dt.timedelta(hours=12), transition_orbit_times[0][-1]),
        (transition_orbit_times[-1][0]-dt.timedelta(hours=12), transition_orbit_times[-1][-1] + dt.timedelta(hours=12)),
                  ]

orbit_data = parse_spice_downsampled([plotting_times[0][0], plotting_times[-1][-1]])

time_series_orbit_data = orbit_data["UTC", "|R|"]

time_series_plot = TimeseriesPanel(time_series_orbit_data)

time_series_plot.plot(show=False)

for i in plotting_times:
    t_start = i[0].isoformat()
    t_end = i[-1].isoformat()
    plots = [PlanarplotPanel(i, plane="X-Y"), PlanarplotPanel(i, plane="X-Z")]
    for plot in plots:
        plot.ax_set_params={
                "title": f"{t_start[10:]} - {t_end[10:]}"
                }
        fig, ax = plot.plot(show=False)
        ax.set_xlim(-4,4)
        ax.set_ylim(-6,4)



delta_t_between_orbits_hist = HistogramPanel(delta_t_numeric, bins=time_bins, minmax=True)

delta_t_between_orbits_hist.ax_set_params = {
        "title": f"Length of all Orbits ({crossing_times[0]} - {crossing_times[-1]})",
        "xlabel": "Orbit length (hours)",
        "ylabel": "Number of orbits",
        "yscale": "log",
        }


delta_t_between_orbits_hist.plot(show=False)
plt.savefig(img_dir + "hollman_orbit_times.svg")
plt.show()
