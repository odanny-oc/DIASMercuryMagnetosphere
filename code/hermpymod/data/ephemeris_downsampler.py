import astropy.units as u
import numpy as np
from astropy.table import QTable, vstack
from hermpy.plotting import Panel
from hermpy.net import ClientSPICE
from astropy.time import Time
import spiceypy as spice
from hermpy.utils import Constants as c
import datetime as dt
import os
import warnings
import subprocess

home_dir = os.getenv("HOME")
data_dir = os.path.join(home_dir, ".ephemeris_data/")
os.makedirs(data_dir, exist_ok = True)

url = "https://zenodo.org/records/17814795/files/hollman_2025_crossing_list.csv?download=1"

save_path = os.path.join(data_dir, "hollman_2025_crossing_list.csv")


EPHEMERIS_FILE = data_dir + 'orbit_ephermis_data_downsampled.ecsv'

# Define custom Mercury radii unit
mercury_rad = c.MERCURY_RADIUS.to("km")
R_M = u.def_unit("R_M", mercury_rad)

u.add_enabled_units(R_M)

mission_start = Time("2011-03-22 00:00:00").to_datetime()
mission_end = Time("2015-04-30 19:00:00").to_datetime()
delta_time = dt.timedelta(hours=12)
sample_number = (mission_end - mission_start) // delta_time


def time_array(start_time, end_time, resolution):
    total_time = (end_time - start_time).total_seconds()
    times = [start_time + dt.timedelta(seconds=s) for s in np.linspace(-1, total_time, int(total_time/resolution))]
    return times


def abs_r(mag_data):
    return np.sqrt(mag_data[0]**2 + mag_data[1]**2 + mag_data[2]**2)


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


def build_ephemeris_table(force_rebuild=False):

    if not os.path.exists(save_path):
        subprocess.run(["wget", "-O", save_path, url])

    if os.path.isfile(EPHEMERIS_FILE) and not force_rebuild:
        print("Ephemeris table already exists. Pass force_rebuild=True to regenerate.")
        return

    time_range = time_array(mission_start, mission_end, 60)

    et = spice.datetime2et(time_range)
    position, _ = spice.spkpos("MESSENGER", et, "MSGR_MSO", "NONE", "Mercury")

    X_MSO = [pos[0] for pos in position]
    Y_MSO = [pos[1] for pos in position]
    Z_MSO = [pos[2] for pos in position]
    
    print('Loaded data')

    # Calculate distance from Mercury
    distance_to_mercury = [abs_r(position[j]) for j in range(len(position))]

    time_range = Time(time_range)
    
    # Create one list of distance from Mercury
    
    print('Calculated |R|')
    
    ephemeris_data = QTable({
                           "UTC": time_range,
                           "|R|": distance_to_mercury * u.Unit("km"),
                           "X MSO": X_MSO * u.Unit("km"),
                           "Y MSO": Y_MSO * u.Unit("km"),
                           "Z MSO": Z_MSO * u.Unit("km"),
                           })
    
    # Save the table
    if force_rebuild:
        ephemeris_data.write(data_dir + 'orbit_ephermis_data_downsampled.ecsv', overwrite = True)
    else:
        ephemeris_data.write(data_dir + 'orbit_ephermis_data_downsampled.ecsv')
    
    print(f'Data Saved')



"""
Fetches, downsampled SPICE data from time range given in requested units. Downloads missing data if needed
"""
def parse_spice_downsampled(time_range = [str, str], units="R_M"):
        t_start, t_end = Time(time_range)
        time = Time(time_range).to_datetime()

        # Handle execptions
        if time[0] < mission_start and time[1] > mission_end:
            raise ValueError(f"Invalid time range given (must lie within {Time(mission_start)} and {Time(mission_end)}).")

        elif time[0] < mission_start:
            time[0] = mission_start
            warnings.warn(f"Start time before mission start ({Time(mission_start)}), starting from mission start")

        elif time[1] > mission_end:
            time[1] = mission_end
            warnings.warn(f"End time after mission end ({Time(mission_end)}), ending at mission end")
        
        if not os.path.isfile(EPHEMERIS_FILE):
            print("Ephemeris table not found — building now...")
            with spice_client.KernelPool():
                build_ephemeris_table()

        table = QTable.read(EPHEMERIS_FILE)

        # O(log n) time lookup via numpy searchsorted
        times_unix = table["UTC"].unix
        i0 = np.searchsorted(times_unix, t_start.unix, side="left")
        i1 = np.searchsorted(times_unix, t_end.unix,   side="right")
        table = table[i0:i1]

        if units != "km":
            for col in table.keys():
                if col == "UTC":
                    continue
                else:
                    table[col] = table[col].to(units)

        return table
