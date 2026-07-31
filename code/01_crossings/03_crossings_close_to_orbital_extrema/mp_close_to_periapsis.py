import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import datetime as dt

from astropy.table import QTable
from astropy.time import Time 

from hermpymod.functions.ephemeris_downsampler import parse_apoapsis_data, parse_periapsis_data, parse_crossing_list
from hermpymod.functions.downsampled_positional_data import parse_spice_downsampled
from hermpymod.functions.encounters import parse_encounters_list
from hermpymod.functions.mag_data_plotter import mag_data_plotter
from hermpymod.functions.plot_all import plot_all_ephemeris
from hermpymod.classes.panels import HistogramPanel


# Load crossing and encounter data set for plots
crossing_data = parse_crossing_list()
encounters_data = parse_encounters_list()


# Load apoapsis and periapsis data
periapsis_data = parse_periapsis_data()
periapsis_times = Time(periapsis_data["UTC"]).to_datetime()

apoapsis_data = parse_apoapsis_data()
apoapsis_times = Time(apoapsis_data["UTC"]).to_datetime()

# Define time range of plot
time_range = [dt.datetime(2013,11,30, 13), dt.datetime(2013,11,30, 21)]

# Plot all planes X-Y... and cylindrical coords of time range
fig, ax = plot_all_ephemeris(time_range, crossings=crossing_data, encounters=encounters_data, add_legend=False)

mag_data_plotter(time_range)

handles, _ = ax[0].get_legend_handles_labels()
fig.legend(handles=handles)

plt.show()
