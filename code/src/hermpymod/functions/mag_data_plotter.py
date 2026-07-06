import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

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

import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from hermpymod.functions.ephemeris_downsampler import parse_crossing_list
from hermpymod.functions.crossings_plotter import plot_crossings

c = ClientMESSENGER()


def total_mag_field(data):
    total_data = data["UTC", "Bx", "By", "Bz"]
    total_mag_data = np.sqrt(total_data["Bx"]**2 + total_data["By"]**2 + total_data["Bz"]**2)
    mag_field_data = QTable(data=[data["UTC"], total_mag_data], names=["UTC", "|B|"])
    return total_data, mag_field_data


def mag_data_plotter(time):
    t_start = time[0]
    t_end = time[-1]

    str_times = [t_start.isoformat(), t_end.isoformat()]
    time_range = TimeRange(t_start, t_end)

    c.query(time_range, "MAG")
    mag_data_encounter = c.fetch()

    mag_table : QTable = parse_messenger_mag(mag_data_encounter, time_range)

    directional_mag_data, total_mag_data = total_mag_field(mag_table)
    totol_mag_plot = TimeseriesPanel(total_mag_data)
    directional_mag_plot = TimeseriesPanel(directional_mag_data)
    
    mag_plot = directional_mag_plot + totol_mag_plot

    fig_mag, ax_mag = mag_plot.plot(show=False)

    fig_mag.suptitle(f"MESSENGER MAG data, taken from {str_times[0][:10]}-{str_times[-1][:10]}")
    lines = ax_mag[1].get_lines()
    lines[0].set_color('k')
    ax_mag[0].axhline(0, ls='--', color='k', label='Zero line')
    ax_mag[0].set_xlabel("Time (UTC)")

    plot_crossings(t_start, t_end, ax_mag)

    for ax in ax_mag:
        ax.legend()

    return fig_mag, ax_mag
