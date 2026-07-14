import astropy.units as u
import numpy as np
from astropy.table import QTable, vstack, hstack
from hermpy.net import ClientSPICE
from astropy.time import Time
import spiceypy as spice
from hermpy.utils import Constants as c # Need for 'Mercury Radii' unit
import datetime as dt
import os
import warnings
import subprocess


home_dir = os.getenv("HOME")
data_dir = os.path.join(home_dir, ".ephemeris_data/")
os.makedirs(data_dir, exist_ok = True)

save_path = os.path.join(data_dir, "hollman_2025_crossing_list.csv")

EPHEMERIS_FILE = data_dir + 'orbit_ephermis_data_downsampled.ecsv'

# 2011-03-23T23:48 first crossing
mission_start = Time("2011-03-18 00:00:00").to_datetime()
mission_end = Time("2015-04-30 19:00:00").to_datetime()

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


def time_array(start_time, end_time, resolution):
    total_time = (end_time - start_time).total_seconds()
    times = [start_time + dt.timedelta(seconds=s) for s in np.linspace(-1, total_time, int(total_time/resolution))]
    return times


def abs_r(mag_data):
    return np.sqrt(mag_data[0]**2 + mag_data[1]**2 + mag_data[2]**2)


"""
Downloads list of peaks with their time and position from Zenodo
"""

def parse_peak_data(force_rebuild=False):

    save_path = os.path.join(data_dir, "peaks_data.csv")

    if force_rebuild:
        os.remove(save_path)

    try:
        peaks_data = QTable.read(save_path)
    except FileNotFoundError:

        url = 'https://zenodo.org/records/21339364/files/peaks_data.csv?download=1'

        subprocess.run(["wget", "-O", save_path, url])

        peaks_data = QTable.read(save_path)

    return peaks_data


"""
Downloads Hollman 2025 crossing list and returns it
"""

def parse_crossing_list(force_rebuild=False):

    crossing_list_dir = os.path.join(data_dir, 'hollman_2025_crossing_list.ecsv')

    save_path = os.path.join(data_dir, "hollman_2025_crossings.csv")


    if force_rebuild:
        os.remove(save_path)
        os.remove(crossing_list_dir)

    # Get crossings if not already downloaded
    try:
        crossing_list = QTable.read(crossing_list_dir)
    except FileNotFoundError:

        url = 'https://zenodo.org/records/20931898/files/hollman_2026_crossing_list.csv?download=1'

        subprocess.run(["wget", "-O", save_path, url])

        crossing_list = QTable.read(save_path)

        time_range = Time(crossing_list["Time"]).to_datetime()
        
        # Get crossing positions
        with spice_client.KernelPool():
            position_table = parse_spice(time_range)

        crossing_table = hstack([position_table, crossing_list["Label"], crossing_list['Trajectory Direction']])
        crossing_table.write(os.path.join(data_dir, "hollman_2025_crossing_list.ecsv"))

        crossing_list = crossing_table


    return crossing_list


"""
Download SPICE data full
"""

def parse_spice(time_range, units="Mercury Radii"):

        time = Time(time_range).to_datetime()
        # Handle execptions
        if time[0] < mission_start and time[-1] > mission_end:
            raise ValueError(f"Invalid time range given (must lie within {Time(mission_start)} and {Time(mission_end)}).")

        # elif time[0] < mission_start:
        #     raise ValueError(f"Start time before mission start ({Time(mission_start)})")

        elif time[-1] > mission_end:
            raise ValueError(f"End time after mission end ({Time(mission_end)})")
        
        # Call data
        et = spice.datetime2et(time_range)
        position, _ = spice.spkpos("MESSENGER", et, "MSGR_MSO", "NONE", "Mercury")

        print('Loaded data')

        X_MSO = [pos[0] for pos in position]
        Y_MSO = [pos[1] for pos in position]
        Z_MSO = [pos[2] for pos in position]

        distance_to_mercury = [abs_r(position[j]) for j in range(len(position))]

        print('Calculated |R|')

        table = QTable({
            "UTC" : Time(time),
            "|R|" : distance_to_mercury * u.Unit("km"),
            "X MSO" : X_MSO * u.Unit("km"),
            "Y MSO" : Y_MSO * u.Unit("km"),
            "Z MSO" : Z_MSO * u.Unit("km"),
            })

        if units != "km":
            for col in table.keys():
                if col == "UTC":
                    continue
                else:
                    table[col] = table[col].to(units)

        return table

"""
Obtains downsampled positional data from SPICE for entire MESSENGER mission and saves it
"""


def build_ephemeris_table(force_rebuild=False):

    _ = parse_crossing_list()

    if os.path.isfile(EPHEMERIS_FILE) and not force_rebuild:
        print("Ephemeris table already exists. Pass force_rebuild=True to regenerate.")
        return

    time_range = time_array(mission_start, mission_end, 60)

    ephemeris_data = parse_spice(time_range, units='km')

    # Save the table
    if force_rebuild:
        ephemeris_data.write(data_dir + 'orbit_ephermis_data_downsampled.ecsv', overwrite = True)
    else:
        ephemeris_data.write(data_dir + 'orbit_ephermis_data_downsampled.ecsv')
    
    print(f'Data Saved')
