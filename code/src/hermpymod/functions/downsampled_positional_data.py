import numpy as np
import os
from astropy.table import QTable
from astropy.time import Time
from hermpy.net import ClientSPICE
from hermpymod.functions.ephemeris_downsampler import build_ephemeris_table
import warnings



home_dir = os.getenv("HOME")
data_dir = os.path.join(home_dir, ".ephemeris_data/")
os.makedirs(data_dir, exist_ok = True)

save_path = os.path.join(data_dir, "hollman_2025_crossing_list.csv")

EPHEMERIS_FILE = data_dir + 'orbit_ephermis_data_downsampled.ecsv'

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


try :
    orbits_downsampled_table = QTable.read(EPHEMERIS_FILE)
except FileNotFoundError:

    with spice_client.KernelPool():
       build_ephemeris_table()

    orbits_downsampled_table = QTable.read(EPHEMERIS_FILE)


"""
Fetches, downsampled SPICE data from time range given in requested units.
"""

def parse_spice_downsampled(time_range = None, units="Mercury Radii"):
        if time_range == None:
            time = [mission_start, mission_end]
            t_start, t_end = Time(time)
        else:
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
        

        table=orbits_downsampled_table

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

