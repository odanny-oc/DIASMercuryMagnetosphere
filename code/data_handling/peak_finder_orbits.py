import numpy as np
import matplotlib.pyplot as plt
from astropy.table import QTable
from sunpy.time import TimeRange
from astropy.table import vstack
from astropy.time import Time

from hermpy.data import parse_messenger_fips, parse_messenger_mag
from hermpy.net import ClientMESSENGER
from hermpy.plotting import MultiPanel, SpectrogramPanel, TimeseriesPanel

import matplotlib.pyplot as plt
import datetime as dt
from datetime import timedelta
from scipy.signal import find_peaks
import spiceypy as spice

import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from hermpymod.data.ephemeris_downsampler import parse_spice_downsampled
from hermpymod.classes.panels import PlanarplotPanel

from hermpy.net import ClientSPICE
from hermpy.utils import Constants as c

import matplotlib

home_dir = os.getenv('HOME')
data_dir = os.path.join(home_dir, ".ephemeris_data/")

spice_client = ClientSPICE()

spice_client.KERNEL_LOCATIONS.update(
    {
        "MESSENGER Frames (tf)": {
            "BASE": "https://naif.jpl.nasa.gov/pub/naif/",
            "DIRECTORY": "pds/data/mess-e_v_h-spice-6-v1.0/messsp_1000/data/fk/",
            "PATTERNS": ["msgr_dyn_v600.tf"],
        },
        "MESSENGER": {
            "BASE": "https://naif.jpl.nasa.gov/pub/naif/",
            "DIRECTORY": "pds/data/mess-e_v_h-spice-6-v1.0/messsp_1000/data/spk/",
            "PATTERNS": ["msgr_??????_??????_??????_od431sc_2.bsp"],
        },
    }
)

mercury_rad = c.MERCURY_RADIUS.to("km").value

matplotlib.use("QtAgg")

ephemeris_data = QTable.read(data_dir + "orbit_ephermis_data_downsampled.ecsv")

crossing_list = QTable.read(data_dir + "hollman_2025_crossing_list.csv")

crossing_times = Time(crossing_list["Time"]).to_datetime()


def abs_r(mag_data):
    return np.sqrt(mag_data[0] ** 2 + mag_data[1] ** 2 + mag_data[2] ** 2)


fig, ax = plt.subplots()

# Plot ephemeris data and create large list of all data
time = ephemeris_data["UTC"].to_datetime()
distance = np.array(ephemeris_data["|R|"])
distance /= mercury_rad
ax.scatter(time, distance, s=0.1, color="C0", rasterized=True)

# Find periapsis to define orbits
peaks = find_peaks(-distance, plateau_size=1,distance=100, height=-1.5)[0]
print(len(peaks))

peaks_data = QTable(
    {
        "UTC": time[peaks],
        "|R|": distance[peaks],
    }
)

average_distance_from_mercury = np.mean(peaks_data["|R|"])

std = np.std(peaks_data["|R|"])

print(average_distance_from_mercury)
print(std)


row_remove_mask = []
for i in range(len(peaks_data) - 1):
    if peaks_data["UTC"][i + 1] - peaks_data["UTC"][i] < timedelta(hours=6):
        row_remove_mask.append(i)
    else:
        continue

for i in row_remove_mask:
    peaks_data.remove_row(i)


# Plot periapsis
ax.scatter(
    time[peaks],
    distance[peaks],
    s=50,
    color="r",
    marker="x",
    label="Peak (Orbit start/end)",
)

time_mask = crossing_times <= time[-1]
crossing_times = crossing_times[time_mask]
crossing_list = crossing_list[time_mask]

with spice_client.KernelPool():
    et = spice.datetime2et(crossing_times)
    crossing_positions, _ = spice.spkpos("MESSENGER", et, "MSGR_MSO", "NONE", "Mercury")

crossing_positions /= mercury_rad

crossing_positions_abs = [
    abs_r(crossing_positions[i]) for i in range(len(crossing_positions))
]

crossing_positions_abs = np.array(crossing_positions_abs)

bs_out_mask = crossing_list["Label"] == "BS_OUT"
bs_in_mask = crossing_list["Label"] == "BS_IN"

mp_out_mask = crossing_list["Label"] == "MP_OUT"
mp_in_mask = crossing_list["Label"] == "MP_IN"

ax.scatter(
    crossing_times[mp_in_mask],
    crossing_positions_abs[mp_in_mask],
    s=50,
    color="b",
    marker="x",
    label="MP_IN",
)
ax.scatter(
    crossing_times[mp_out_mask],
    crossing_positions_abs[mp_out_mask],
    s=50,
    color="purple",
    marker="x",
    label="MP_OUT",
)

ax.scatter(
    crossing_times[bs_in_mask],
    crossing_positions_abs[bs_in_mask],
    s=50,
    color="orange",
    marker="x",
    label="BS_IN",
)
ax.scatter(
    crossing_times[bs_out_mask],
    crossing_positions_abs[bs_out_mask],
    s=50,
    color="yellow",
    marker="x",
    label="BS_OUT",
)

ax.set_title(f"MESSENGER Distance from Mercury")
ax.set_ylabel(r"Distance from Mercury |R| (Mercury Radii $R_M$)")
ax.set_xlabel(f"Time YY:MM:DD (UTC)")
ax.grid()

plt.savefig(data_dir + "orbits_plot_with_peaks.svg")

# Save peak data
peaks_data.write(data_dir + "peaks_data.csv", overwrite=True)

plt.legend()
print("Plot finished!")

# 2D Planar Plots

plane_plot = PlanarplotPanel(["2012-03-15", "2012-03-18"], plane="X-Y")
fig, ax = plane_plot.plot(show=False)

plane_time_total = plane_plot.time
plane_time = [plane_time_total[0], plane_time_total[-1]]

ax.set_xlabel(r"X MSO (Mercury Radii $R_M$)")
ax.set_ylabel(r"Y MSO (Mercury Radii $R_M$)")

ax.set_title("MESSENGER Trajectory X-Y Plane (MSO)")

def match_times(time_array, compare_array):
    mask = (time_array[0] <= compare_array) & ((time_array[-1] >= compare_array))
    return mask

plane_mask = match_times(plane_time, crossing_times)

crossing_list = crossing_list[plane_mask]
crossing_positions = crossing_positions[plane_mask]

bs_out_mask = crossing_list["Label"] == "BS_OUT"
bs_in_mask = crossing_list["Label"] == "BS_IN"

mp_out_mask = crossing_list["Label"] == "MP_OUT"
mp_in_mask = crossing_list["Label"] == "MP_IN"

ax.scatter(
    crossing_positions[:,0][mp_in_mask],
    crossing_positions[:,1][mp_in_mask],
    s=50,
    color="b",
    marker="x",
    label="MP_IN",
)
ax.scatter(
    crossing_positions[:,0][mp_out_mask],
    crossing_positions[:,1][mp_out_mask],
    s=50,
    color="purple",
    marker="x",
    label="MP_OUT",
)

ax.scatter(
    crossing_positions[:,0][bs_in_mask],
    crossing_positions[:,1][bs_in_mask],
    s=50,
    color="orange",
    marker="x",
    label="BS_IN",
)
ax.scatter(
    crossing_positions[:,0][bs_out_mask],
    crossing_positions[:,1][bs_out_mask],
    s=50,
    color="yellow",
    marker="x",
    label="BS_OUT",
)

ax.set_aspect('equal')
plt.show()
