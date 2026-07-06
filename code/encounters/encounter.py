import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from astropy.table import QTable, vstack, Column
import astropy.units as u
from astropy.time import Time
from astropy.table import QTable, hstack, vstack
from sunpy.time import TimeRange

from hermpy.plotting import Panel, TimeseriesPanel
from hermpy.data import parse_messenger_fips, parse_messenger_mag
from hermpy.net import ClientSPICE, ClientMESSENGER
from hermpy.utils import Constants as c

import spiceypy as spice

import datetime as dt

import json
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from hermpymod.classes.panels import PlanarplotPanel, HistogramPanel
from hermpymod.functions.ephemeris_downsampler import parse_spice_downsampled, parse_crossing_list
from hermpymod.functions.crossings_plotter import plot_crossings
from hermpymod.functions.mag_data_plotter import mag_data_plotter


c = ClientMESSENGER()

home_dir = os.getenv('HOME')
data_dir = os.path.join(home_dir, '.ephemeris_data/')

plt.style.use("../../plots_and_images/presentation.mplstyle")

peaks_data = QTable.read(data_dir + 'peaks_data.csv')

crossing_list = parse_crossing_list()

crossing_label = crossing_list["Label"]
crossing_direction = crossing_list["Trajectory Direction"]
crossing_time = crossing_list["UTC"].to_datetime()

encounters_data = []
time_start = []
time_end = []
label = []
cross_number = []
dts= []
encounter_dt = []


for i in range(len(crossing_time) - 1):
    time_delta = crossing_time[i+1] - crossing_time[i]
    crossing_type = crossing_label[i][:2]
    if crossing_type in crossing_label[i+1] and crossing_direction[i] == crossing_direction[i+1] and time_delta<=dt.timedelta(hours=13):
        cross_number.append(crossing_list[i])
        dts.append(time_delta)
    else:
        cross_number.append(crossing_list[i])
        dts.append(time_delta)
        time_start.append(cross_number[0]["UTC"].iso)
        time_end.append(cross_number[-1]["UTC"].iso)
        label.append(crossing_type + crossing_direction[i][0])
        rows = vstack(cross_number)

        encounters_data.append(rows)
        encounter_dt.append(dts)
        cross_number = []
        dts =[]

encounter_duration = [(Time(time_end[i]).to_datetime() - Time(time_start[i]).to_datetime()).total_seconds()/3600 for i in range(len(time_start))]

encounters_list = {
        "Time Start": time_start,
        "Time End": time_end,
        "Label": label,
        "Encounter Duration": encounter_duration
        }

print("Number of Encounters", len(encounters_list["Time Start"]))
hollman_encounters_list = pd.DataFrame(encounters_list)

times_dt = Time(encounters_list["Time Start"]).to_datetime()
i0 = np.searchsorted(times_dt, dt.datetime(2012,10,14, 6), side="left")


hollman_encounters_list.to_csv(data_dir + "hollman_encounters_list_2025.csv", index=False)


for i in range(len(encounters_data)):
    encounters_data[i].add_column(Column(encounter_dt[i], name='dt'))


def total_mag_field(data):
    total_data = data["UTC", "Bx", "By", "Bz"]
    total_mag_data = np.sqrt(total_data["Bx"]**2 + total_data["By"]**2 + total_data["Bz"]**2)
    mag_field_data = QTable(data=[data["UTC"], total_mag_data], names=["UTC", "|B|"])
    return total_data, mag_field_data


panels = []
mag_encounters = [i for i in encounters_data if len(i) >= 2]
mag_encounters = [i for i in mag_encounters if len(i) %2 == 1]

print(len(mag_encounters))
slicing= len(mag_encounters)//3

mag_encounters = mag_encounters[::slicing]
mag_encounters.append(encounters_data[i0])

# MAG plots for encounters

test_times = []

for i in mag_encounters:
    encounter_time = i["UTC"].to_datetime()
    timewindow = dt.timedelta(minutes=300)
    t_start, t_end = encounter_time[0] - timewindow, encounter_time[-1] + timewindow

    # Plot MAG data in time interval and all crossings
    fig_mag, ax_mag = mag_data_plotter([t_start, t_end])

    # Plot encounter on graph
    ax_mag[0].axvspan(encounter_time[0], encounter_time[-1], alpha=0.5, color='orange', label='Encounter')
    ax_mag[1].axvspan(encounter_time[0], encounter_time[-1], alpha=0.5, color='orange', label='Encounter')

    plane_plot_xy = PlanarplotPanel([t_start, t_end], plane='X-Y', crossings=True)
    plane_plot_yz = PlanarplotPanel([t_start, t_end], plane='Y-Z', crossings=True)
    plane_plot =[ plane_plot_xy , plane_plot_yz]

    # Spatial Plots
    for plane in plane_plot:
        _, ax_pos = plane.plot(show=False)
        encounter_mask = (plane.time_range >= encounter_time[0]) & (plane.time_range <= encounter_time[-1])
        plane_positions = plane.positions()[encounter_mask]
        ax_pos.plot(plane_positions[plane.labels[0]], plane_positions[plane.labels[1]], zorder =2, color ='orange', lw=3)
        ax_pos.set_title(f"MESSENGER SPICE data, taken from {t_start}-{t_end}")
        ax_pos.set_aspect('equal')
        ax_pos.legend()

plt.show()

# with open("test_times", "w") as f:
#     json.dump(test_times, f)
print(test_times)

mp_encounters = [i for i in encounters_data if "MP" in i[0]['Label']]
bs_encounters = [i for i in encounters_data if "BS" in i[0]['Label']]

number_of_crossing_encounter = [len(i) for i in encounters_data]
number_of_crossing_encounter_mp = [len(i) for i in mp_encounters]
number_of_crossing_encounter_bs = [len(i) for i in bs_encounters]

encounter_bins = np.arange(0, 17, 1)
encounters_hist = HistogramPanel(number_of_crossing_encounter, bins=encounter_bins, color='red')

encounters_hist.ax_set_params = {
        "title": "Number of crossings per 'encounter' (defined by crossings of the same type that happen within ~20mins)",
        "xlabel": "Number of crossings",
        "ylabel": "Number of encounters",
        "yscale": "log"
        }

encounters_hist_mp = HistogramPanel(number_of_crossing_encounter_mp, bins=encounter_bins, color='purple')

encounters_hist_mp.ax_set_params = {
        "title": "Number of MP crossings per 'encounter' (defined by crossings of the same type that happen within ~20mins)",
        "xlabel": "Number of crossings",
        "ylabel": "Number of encounters",
        "yscale": "log"
        }

encounters_hist_bs = HistogramPanel(number_of_crossing_encounter_bs, bins=encounter_bins, color='yellow')

encounters_hist_bs.ax_set_params = {
        "title": "Number of BS crossings per 'encounter' (defined by crossings of the same type that happen within ~20mins)",
        "xlabel": "Number of crossings",
        "ylabel": "Number of encounters",
        "yscale": "log"
        }

encounters_sub_hist = encounters_hist_bs + encounters_hist_mp

encounters_hist.plot(show=False)
encounters_sub_hist.plot(show=False)

plt.show()
