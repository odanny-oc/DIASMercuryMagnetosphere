import numpy as np
import matplotlib.pyplot as plt
from astropy.table import QTable
from sunpy.time import TimeRange
from astropy.time import Time

import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import spiceypy as spice

import os

from hermpymod.classes.panels import PlanarplotPanel
from hermpymod.functions.ephemeris_downsampler import parse_crossing_list, parse_spice_downsampled

from hermpy.net import ClientSPICE
from hermpy.utils import Constants as c

import matplotlib as mpl

home_dir = os.getenv('HOME')
data_dir = os.path.join(home_dir, ".ephemeris_data/")
images_dir = "../../plots_and_images/"

mpl.use("QtAgg")

plt.style.use(images_dir + "presentation.mplstyle")

ephemeris_data = parse_spice_downsampled()

crossing_list = parse_crossing_list()

crossing_times = Time(crossing_list["UTC"]).to_datetime()


def abs_r(mag_data):
    return np.sqrt(mag_data[0] ** 2 + mag_data[1] ** 2 + mag_data[2] ** 2)


fig, ax = plt.subplots()

# Plot ephemeris data and create large list of all data
time = ephemeris_data["UTC"].to_datetime()
distance = np.array(ephemeris_data["|R|"])

ax.plot(time, distance, color="C0",lw =0.4, zorder=1)

# Find periapsis to define orbits
peaks = find_peaks(-distance, plateau_size=1,distance=100, height=-1.5)[0]

# Include last partial orbit
peaks = np.append(peaks, -1)
# peaks = np.insert(peaks, 0 , 0)

print("Number of peaks", len(peaks))

delta_t_between_orbits = [0]

peak_times = Time(time[peaks]).to_datetime()

orbit_number = [0] 


for i in range(len(peak_times) - 1):
    time_start =  peak_times[i]
    time_end = peak_times[i + 1]
    orbit_number.append(i+1)
    orbit_time = time_end - time_start 
    delta_t_between_orbits.append(orbit_time)


peaks_data = QTable(
    {
        "UTC": time[peaks],
        "X MSO": ephemeris_data["X MSO"][peaks],
        "Y MSO": ephemeris_data["Y MSO"][peaks],
        "Z MSO": ephemeris_data["Z MSO"][peaks],
        "|R|": distance[peaks],
        "Delta t": delta_t_between_orbits,
        "Orbit Number": orbit_number,
    }
)

# Plot periapsis
ax.scatter(
    time[peaks],
    distance[peaks],
    s=50,
    color="r",
    marker="x",
    label="Peak (Orbit start/end)",
)

crossing_positions = crossing_list["X MSO", "Y MSO", "Z MSO"]

crossing_positions_abs = [
    abs_r(crossing_positions[i]) for i in range(len(crossing_positions))
]

# Convert to numpy array to use mask
crossing_positions_abs = np.array([i.value for i in crossing_positions_abs])

crossings_dict ={
        "BS_OUT": {"mask": crossing_list["Label"] == "BS_OUT", "color": "yellow"},
        "BS_IN": {"mask": crossing_list["Label"] == "BS_IN", "color": "orange"},
        "MP_OUT": {"mask": crossing_list["Label"] == "MP_OUT", "color": "purple"},
        "MP_IN": {"mask": crossing_list["Label"] == "MP_IN", "color": "blue"},
        }


for label in crossings_dict.keys():
    dict = crossings_dict[label]
    mask = dict["mask"]
    color = dict["color"]
    ax.scatter(
        crossing_times[mask],
        crossing_positions_abs[mask],
        s=50,
        color=color,
        marker="x",
        label=label,
    )

ax.set_title(f"MESSENGER Distance from Mercury ({time[peaks][0]} - {time[peaks][-1]})")
ax.set_ylabel(r"Distance from Mercury $\|R\|$ (Mercury Radii $R_M$)")
ax.set_xlabel(f"Time MM:DD:HH (UTC)")
# ax.set_xlim(dt.datetime(2011,3,31), dt.datetime(2011,4,3))
ax.grid()

plt.legend(loc='upper left')

# plt.savefig(images_dir + "orbits_plot_with_peaks.svg")

# Save peak data
peaks_data.write(data_dir + "peaks_data.csv", overwrite=True)

print("Plot finished!")

# 2D Planar Plots

times = ["2011-03-30", "2011-04-30"]
planes = ["X-Y", "X-Z", "Y-Z"]

mpl.rcParams['path.simplify'] = True
mpl.rcParams['path.simplify_threshold'] = 1.0

for plane in planes:
    plane_plot_mp = PlanarplotPanel(times, plane=plane, crossings=True, BS=False, alpha=0.05)
    plane_plot_bs = PlanarplotPanel(times, plane=plane, crossings=True, MP=False, alpha=0.05)

    plane_plot = plane_plot_bs + plane_plot_mp

    fig, ax = plane_plot.plot(show=False)
    fig.suptitle(f"MESSENGER Trajectory, {times}, {plane} Plane")
    # plt.savefig(images_dir + f"orbit_{plane}_crossings.svg")

plt.show()
