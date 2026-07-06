import numpy as np
from astropy.table import QTable
from sunpy.time import TimeRange
from astropy.time import Time

from hermpy.data import parse_messenger_fips, parse_messenger_mag
from hermpy.net import ClientMESSENGER

from datetime import timedelta
import os
import sys

from hermpymod.functions.ephemeris_downsampler import parse_crossing_list

home_dir = os.getenv('HOME')
data_dir = os.path.join(home_dir, ".ephemeris_data/")
img_dir = "../../plots_and_images/"



peak_data = QTable.read(data_dir + "peaks_data.csv")

peak_times = Time(peak_data["UTC"]).to_datetime()

crossing_data = parse_crossing_list()

crossing_numbers = []

crossing_times =Time(crossing_data["UTC"]).to_datetime()

orbit_list = []

delta_t_between_orbits = []


def orbit_data():
    for i in range(len(peak_times) - 1):
        time_start =  peak_times[i]
        time_end = peak_times[i + 1]
        orbit_time = time_end - time_start 
        delta_t_between_orbits.append(orbit_time)
        if orbit_time > timedelta(hours=13) or orbit_time < timedelta(hours=7):
            continue
        else:
            mask = (crossing_times >= peak_times[i]) & (crossing_times <= peak_times[i + 1])
            orbit_list.append(crossing_data[mask])

    return orbit_list
