import os
from pathlib import Path
os.chdir(Path(__file__).resolve().parent)

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import QTable
from sunpy.time import TimeRange
from astropy.time import Time

import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import spiceypy as spice
import datetime as dt

import os

from hermpymod.classes.panels import PlanarplotPanel
from hermpymod.functions.ephemeris_downsampler import parse_crossing_list
from hermpymod.functions.downsampled_positional_data import parse_spice_downsampled

from hermpy.net import ClientSPICE
from hermpy.utils import Constants as c

import matplotlib as mpl

home_dir = os.getenv('HOME')
from hermpymod.paths import DATA_DIR
data_dir = DATA_DIR
images_dir = "../../plots_and_images/"

mpl.use("QtAgg")

plt.style.use(images_dir + "presentation.mplstyle")


def abs_r(mag_data):
    return np.sqrt(mag_data[0] ** 2 + mag_data[1] ** 2 + mag_data[2] ** 2)


ephemeris_data = parse_spice_downsampled()

crossing_list = parse_crossing_list(force_rebuild=True)

crossing_times = Time(crossing_list["UTC"]).to_datetime()

fig, ax = plt.subplots()

# Plot ephemeris data and create large list of all data
time = ephemeris_data["UTC"].to_datetime()
distance = np.array(ephemeris_data["|R|"])

ax.plot(time, distance, color="C0",lw =0.4, zorder=1)

# Find periapsis to define orbits
periapsides = find_peaks(-distance, plateau_size=1,distance=100, height=-1.5)[0]
apoapsides = find_peaks(distance, plateau_size=1,distance=100, height=-1.5)[0]

# Include last partial orbit
periapsides = np.append(periapsides, -1)
# peaks = np.insert(peaks, 0 , 0)

print("Number of periapsides", len(periapsides))
print("Number of apoapsides", len(apoapsides))

delta_t_between_periapsides = [0]

periapsides_times = Time(time[periapsides]).to_datetime()
print(periapsides_times[0])

periapsis_orbit_number = [0] 

sig_figs = 4

for i in range(len(periapsides_times) - 1):
    time_start =  periapsides_times[i]
    time_end = periapsides_times[i + 1]
    periapsis_orbit_number.append(i+1)
    orbit_time = round((time_end - time_start).total_seconds()/3600, sig_figs)
    delta_t_between_periapsides.append(orbit_time)

delta_t_between_apoapsides = [0]

apoapsides_times = Time(time[apoapsides]).to_datetime()

apoapsis_orbit_number = [0] 


for i in range(len(apoapsides_times) - 1):
    time_start =  apoapsides_times[i]
    time_end = apoapsides_times[i + 1]
    apoapsis_orbit_number.append(i+1)
    orbit_time = round((time_end - time_start).total_seconds()/3600, sig_figs)
    delta_t_between_apoapsides.append(orbit_time)


def round_datetime_to_second(t):
    return t.replace(microsecond=0) + dt.timedelta(seconds=round(t.microsecond / 1e6))


periapsis_utc_rounded = [round_datetime_to_second(t) for t in time[periapsides]]
apoapsis_utc_rounded = [round_datetime_to_second(t) for t in time[apoapsides]]


periapsis_data = QTable(
    {
        "UTC": periapsis_utc_rounded,
        "X MSO": round(ephemeris_data["X MSO"][periapsides], sig_figs),
        "Y MSO": round(ephemeris_data["Y MSO"][periapsides], sig_figs),
        "Z MSO": round(ephemeris_data["Z MSO"][periapsides], sig_figs),
        "|R|": np.round(distance[periapsides], sig_figs),
        "Orbit Length": delta_t_between_periapsides,
        "Orbit Number": periapsis_orbit_number,
    }
)


apoapsis_data = QTable(
    {
        "UTC": apoapsis_utc_rounded,
        "X MSO": round(ephemeris_data["X MSO"][apoapsides], sig_figs),
        "Y MSO": round(ephemeris_data["Y MSO"][apoapsides], sig_figs),
        "Z MSO": round(ephemeris_data["Z MSO"][apoapsides], sig_figs),
        "|R|": np.round(distance[apoapsides], sig_figs),
        "Orbit Length": delta_t_between_apoapsides,
        "Orbit Number": apoapsis_orbit_number,
    }
)


# Plot periapsis
ax.scatter(
    time[periapsides],
    distance[periapsides],
    s=50,
    color="r",
    marker="x",
    label="Periapsis (Orbit start/end)",
)


# Plot periapsis
ax.scatter(
    time[apoapsides],
    distance[apoapsides],
    s=50,
    color="green",
    marker="x",
    label="Apoapsis (Orbit start/end)",
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

ax.set_title(f"MESSENGER Distance from Mercury ({time[periapsides][0]} - {time[periapsides][-1]})")
ax.set_ylabel(r"Distance from Mercury $\|R\|$ (Mercury Radii $R_M$)")
ax.set_xlabel(f"Time MM:DD:HH (UTC)")
# ax.set_xlim(dt.datetime(2011,3,31), dt.datetime(2011,4,3))
ax.grid()

plt.legend(loc='upper left')

# plt.savefig(images_dir + "orbits_plot_with_peaks.svg")

# Save peak data
periapsis_data.write(data_dir + "periapsis_data.csv", overwrite=True)
apoapsis_data.write(data_dir + "apoapsis_data.csv", overwrite=True)

print("Plot finished!")

# 2D Planar Plots

times = ["2011-03-30", "2011-04-30"]
planes = ["X-Y", "X-Z", "Y-Z"]

mpl.rcParams['path.simplify'] = True
mpl.rcParams['path.simplify_threshold'] = 1.0

for plane in planes:
    plane_plot_mp = PlanarplotPanel(times, plane=plane, crossings=crossing_list, BS=False, alpha=0.05)
    plane_plot_bs = PlanarplotPanel(times, plane=plane, crossings=crossing_list, MP=False, alpha=0.05)

    plane_plot = plane_plot_bs + plane_plot_mp

    fig, ax = plane_plot.plot(show=False)
    fig.suptitle(f"MESSENGER Trajectory, {times}, {plane} Plane")
    # plt.savefig(images_dir + f"orbit_{plane}_crossings.svg")

plt.show()
