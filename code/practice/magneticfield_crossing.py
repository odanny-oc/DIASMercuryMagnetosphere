import numpy as np
import matplotlib.pyplot as plt
from astropy.table import QTable
from sunpy.time import TimeRange
from astropy.time import Time
import random

from hermpy.data import parse_messenger_fips, parse_messenger_mag
from hermpy.net import ClientMESSENGER
from hermpy.plotting import MultiPanel, SpectrogramPanel, TimeseriesPanel

import matplotlib.pyplot as plt
from matplotlib import colormaps as cm
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib.colors import Normalize
from PIL import Image


"""
Function to calculate the total magnetic field strenght as well as extract only the magnetic field data from the MESSENGER data set.
"""
def total_mag_field(data):
    total_data = data["UTC", "Bx", "By", "Bz"]
    total_mag_data = np.sqrt(total_data["Bx"]**2 + total_data["Bx"]**2 + total_data["Bz"]**2)
    mag_field_data = QTable(data=[data["UTC"], total_mag_data], names=["UTC", "|B|"])
    return total_data, mag_field_data


c = ClientMESSENGER()
time_range_1_hour = TimeRange("2011-09-27T12:00", "2011-09-27T13:00")
time_range_2_hour = TimeRange("2011-09-27T12:00", "2011-09-27T14:00")
time_range_8_hour = TimeRange("2011-09-27T12:00", "2011-09-27T20:00")

mag_data = {}

c.query(time_range_1_hour, "MAG")
all_mag_data_1_hour = c.fetch()

c.query(time_range_2_hour, "MAG")
all_mag_data_2_hour = c.fetch()

c.query(time_range_8_hour, "MAG")
all_mag_data_8_hour = c.fetch()

# Dictionary of each time interval
mag_data[0]: QTable = parse_messenger_mag(all_mag_data_1_hour, time_range_1_hour)
mag_data[1]: QTable = parse_messenger_mag(all_mag_data_2_hour, time_range_2_hour)
mag_data[2]: QTable = parse_messenger_mag(all_mag_data_8_hour, time_range_8_hour)

mag_plotting_data = {}
mag_panels = {}
for i in mag_data:
    mag_plotting_data[i] = total_mag_field(mag_data[i])
    mag_panels[i] = TimeseriesPanel(mag_plotting_data[i][0]) + TimeseriesPanel(mag_plotting_data[i][1])


model_crossing_list = QTable.read(r"../data/hollman_2025_crossing_list.csv")

fig, ax = mag_panels[2].plot(show=False)

xmin = mag_plotting_data[2][0]["UTC"][0].to_datetime()
xmax = mag_plotting_data[2][0]["UTC"][-1].to_datetime()

# Change total magnetic field to black
lines = ax[1].get_lines()
lines[0].set_color('k')

# Find all crossings in our given time interval
data_in_timerange_bool= ["2011-09-27" in t for t in model_crossing_list["Time"]]

true_indices = [i for i in range(len(data_in_timerange_bool)) if data_in_timerange_bool[i]==True]

crossings = model_crossing_list[true_indices]

crossings["Time"] = Time(crossings["Time"])


for time,label in zip(crossings["Time"], crossings["Label"]):
    if "MP" in label:
        if "OUT" in label:
            ax[1].axvline(time.to_datetime(), color = 'purple', label=label)
        else:
            ax[1].axvline(time.to_datetime(), color = 'blue', label=label)
    elif "BS" in label:
        if "OUT" in label:
            ax[1].axvline(time.to_datetime(), color = 'yellow', label=label)
        else:
            ax[1].axvline(time.to_datetime(), color = 'red', label=label)

plt.xlim(xmin, xmax)

handles, labels = plt.gca().get_legend_handles_labels()

by_label = dict(zip(labels, handles))

ax[1].set_yscale('log')

ax[1].legend(by_label.values(), by_label.keys())

plt.show()
