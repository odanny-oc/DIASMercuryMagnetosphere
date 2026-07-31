import numpy as np
from astropy.table import QTable
from sunpy.time import TimeRange
from astropy.time import Time


from datetime import timedelta
import os
from hermpymod.functions.ephemeris_downsampler import parse_crossing_list, parse_periapsis_data


home_dir = os.getenv('HOME')
data_dir = os.path.join(home_dir, ".ephemeris_data/")

peak_data = parse_periapsis_data()

crossing_data = parse_crossing_list()

mask = ["UNPHYSICAL" not in label for label in crossing_data["Label"]]
crossing_data = crossing_data[mask]


def strip_leading_empty(lst):
    i = 0
    while i < len(lst) and not lst[i]:
        i += 1
    return lst[i:]


"""
Returns crossings in lists corresponding to their orbit. Strips leading empty lists, i.e for Hollman list, first list is the 12th orbit
"""

def orbit_data(crossing_data=crossing_data, peak_data=peak_data):
    orbit_list = []
    crossing_times =Time(crossing_data["UTC"]).to_datetime()
    peak_times = Time(peak_data["UTC"]).to_datetime()
    for i in range(len(peak_times) - 1):
        time_start =  peak_times[i]
        time_end = peak_times[i + 1]
        orbit_time = time_end - time_start 
        if orbit_time > timedelta(hours=13) or orbit_time < timedelta(hours=7):
            continue
        else:
            mask = (crossing_times >= peak_times[i]) & (crossing_times <= peak_times[i + 1])
            orbit_number = np.full(len(crossing_data[mask]), peak_data["Orbit Number"][i + 1])
            table = crossing_data[mask].copy()
            table["Orbit Number"] = orbit_number
            orbit_list.append(table)

    orbit_list_stripped = strip_leading_empty(orbit_list)

    return orbit_list_stripped
