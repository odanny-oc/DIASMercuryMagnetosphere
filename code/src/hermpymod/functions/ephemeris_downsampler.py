import numpy as np

import astropy.units as u
from astropy.table import QTable, vstack, hstack
from astropy.time import Time
from hermpy.utils import Constants as c # Need for 'Mercury Radii' unit
from hermpy.net import ClientSPICE

import spiceypy as spice
import datetime as dt
import os
import warnings
import subprocess
from hermpymod.paths import DATA_DIR


Zd = c.DIPOLE_OFFSET.to("Mercury Radii")

data_dir = DATA_DIR

EPHEMERIS_FILE = os.path.join(data_dir , 'orbit_ephermis_data_downsampled.ecsv')

# 2011-03-23T23:48 first crossing
mission_start = Time("2011-03-18 00:50:00").to_datetime()
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



def parse_periapsis_data(force_rebuild=False):
    """
    Downloads list of periapsides with their time and position from Zenodo.

    Caches the result locally as a CSV; subsequent calls read from disk
    instead of re-downloading, unless force_rebuild=True.

    force_rebuild: if True, deletes any existing cached CSV first, forcing
        a fresh download from Zenodo.

    Returns:
        periapsis_data: QTable of periapsis times, positions, orbit durations, and orbit numbers.
    """
    save_path = os.path.join(data_dir, "periapsis_data.csv")

    if force_rebuild:
        try:
            os.remove(save_path)
        except FileNotFoundError:
            pass

    # Load from cache if it exists, otherwise download from Zenodo
    try:
        periapsis_data = QTable.read(save_path)
    except FileNotFoundError:
        url = 'https://zenodo.org/records/21393912/files/periapsis_data.csv?download=1'
        subprocess.run(["curl", "-Lo", save_path, url])
        periapsis_data = QTable.read(save_path)

    return periapsis_data


def parse_apoapsis_data(force_rebuild=False):
    """
    Downloads list of apoapsides with their time and position from Zenodo.

    Caches the result locally as a CSV; subsequent calls read from disk
    instead of re-downloading, unless force_rebuild=True.

    force_rebuild: if True, deletes any existing cached CSV first, forcing
        a fresh download from Zenodo.

    Returns:
        apoapsis_data: QTable of apoapsis times and positions, orbit durations, and orbit numbers.
    """
    save_path = os.path.join(data_dir, "apoapsis_data.csv")

    if force_rebuild:
        try:
            os.remove(save_path)
        except FileNotFoundError:
            pass

    # Load from cache if it exists, otherwise download from Zenodo
    try:
        apoapsis_data = QTable.read(save_path)
    except FileNotFoundError:
        url = 'https://zenodo.org/records/21393912/files/apoapsis_data.csv?download=1'
        subprocess.run(["curl", "-Lo", save_path, url])
        apoapsis_data = QTable.read(save_path)

    return apoapsis_data


def parse_crossing_list(force_rebuild=False):
    """
    Downloads the Hollman 2026 magnetopause/bow shock crossing list and
    returns it, with corresponding spacecraft position appended.

    Two files are cached: the raw downloaded crossing list (save_path) and
    the processed table with SPICE-derived positions merged in
    (crossing_list_dir). Only the processed table is checked/read on
    subsequent calls.

    force_rebuild: if True, deletes both cached files first, forcing a
        fresh download and reprocessing.

    Returns:
        crossing_list: QTable with columns for spacecraft position (from
            SPICE), "Label" (crossing type), and "Trajectory Direction".
    """
    crossing_list_dir = os.path.join(data_dir, 'hollman_2026_crossing_list.csv')
    save_path = os.path.join(data_dir, "hollman_2026_crossings.csv")

    if force_rebuild:
        try:
            os.remove(save_path)
            os.remove(crossing_list_dir)
        except FileNotFoundError:
            pass

    # Get crossings if not already downloaded/processed
    try:
        crossing_list = QTable.read(crossing_list_dir)
    except FileNotFoundError:
        try:
            crossing_list = QTable.read(save_path)
        except FileNotFoundError:
            # Download the raw crossing list from Zenodo
            url = 'https://zenodo.org/records/21392216/files/hollman_2026_crossing_list.csv?download=1'
            subprocess.run(["curl", "-Lo", save_path, url])
            crossing_list = QTable.read(save_path)

        time_range = Time(crossing_list["Time"]).to_datetime()

        # Look up spacecraft position at each crossing time via SPICE
        with spice_client.KernelPool():
            position_table = parse_spice(time_range)

        # Combine position data with crossing label/direction, then cache to disk
        crossing_table = hstack([position_table, crossing_list["Label"], crossing_list['Trajectory Direction']])
        crossing_table.write(crossing_list_dir)

        crossing_list = crossing_table

    return crossing_list


def parse_spice(time_array, units="Mercury Radii", frame="MSO"):
    """
    Queries SPICE for MESSENGER's position relative to Mercury at each given
    time, and returns it as a table in the requested reference frame/units.

    time_array: array-like of times (anything astropy.time.Time can parse)
        to query positions for. Must lie within the MESSENGER mission
        timespan (mission_start to mission_end).
    units: unit to convert the position columns to (e.g. "Mercury Radii",
        "km"). "UTC" column is left untouched.
    frame: coordinate frame to return positions in.
        "MSO" - Mercury Solar Orbital
        "MSM" - Mercury Solar Magnetospheric (MSO shifted along Z by the
                dipole offset Zd)
        "All" - both MSO and MSM columns included

    Returns:
        table: QTable with "UTC", "|R|" (distance to Mercury), and position
            columns for the requested frame(s).
    """
    time = Time(time_array).to_datetime()

    if isinstance(time, np.ndarray):
        # Handle exceptions: requested range must fall within mission dates
        if time[0] < mission_start and time[-1] > mission_end:
            raise ValueError(f"Invalid time range given (must lie within {Time(mission_start)} and {Time(mission_end)}).")

        # elif time[0] < mission_start:
        #     raise ValueError(f"Start time before mission start ({Time(mission_start)})")

        elif time[-1] > mission_end:
            raise ValueError(f"End time after mission end ({Time(mission_end)})")

    else:
        print("arrayed time")

    # Convert to ephemeris time and query spacecraft position relative to Mercury
    et = spice.datetime2et(time)
    position, _ = spice.spkpos("MESSENGER", et, "MSGR_MSO", "NONE", "Mercury")

    print('Loaded data')

    X_MSO = np.array([pos[0] for pos in position])
    Y_MSO = np.array([pos[1] for pos in position])
    Z_MSO = np.array([pos[2] for pos in position])

    distance_to_mercury = [abs_r(position[j]) for j in range(len(position))]

    print('Calculated |R|')

    if frame == "MSO":
        table = QTable({
            "UTC" : Time(time),
            "|R|" : distance_to_mercury * u.Unit("km"),
            "X MSO" : X_MSO * u.Unit("km"),
            "Y MSO" : Y_MSO * u.Unit("km"),
            "Z MSO" : Z_MSO * u.Unit("km"),
            })

    elif frame == "MSM":
        # Shift Z by the dipole offset to convert MSO -> MSM
        distance_to_mercury = abs_r([X_MSO, Y_MSO, Z_MSO + Zd.to("km").value])
        table = QTable({
            "UTC" : Time(time),
            "|R|" : distance_to_mercury * u.Unit("km"),
            "X MSM" : X_MSO * u.Unit("km"),
            "Y MSM" : Y_MSO * u.Unit("km"),
            "Z MSM" : Z_MSO * u.Unit("km") - Zd.to("km"),
            })

    elif frame == "All":
        table = QTable({
            "UTC" : Time(time),
            "|R|" : distance_to_mercury * u.Unit("km"),
            "X MSO" : X_MSO * u.Unit("km"),
            "Y MSO" : Y_MSO * u.Unit("km"),
            "Z MSO" : Z_MSO * u.Unit("km"),
            "X MSM" : X_MSO * u.Unit("km"),
            "Y MSM" : Y_MSO * u.Unit("km"),
            "Z MSM" : Z_MSO * u.Unit("km") - Zd.to("km"),
            })
    else:
        raise ValueError("Invalid frame, not one of 'MSO', 'MSM', or 'All'")

    # Convert all non-UTC columns to the requested output units
    if units != "km":
        for col in table.keys():
            if col == "UTC":
                continue
            else:
                table[col] = table[col].to(units)

    return table


def build_ephemeris_table(force_rebuild=False, frame="All"):
    """
    Builds a downsampled positional ephemeris table (60s cadence, by
    default) for the entire MESSENGER mission from SPICE, and saves it to
    disk at EPHEMERIS_FILE.

    Skips rebuilding if EPHEMERIS_FILE already exists, unless
    force_rebuild=True.

    force_rebuild: if True, regenerates and overwrites the ephemeris table
        even if one already exists on disk.
    frame: coordinate frame to pass through to parse_spice ("MSO", "MSM",
        or "All").

    Returns:
        None. Writes the ephemeris table to EPHEMERIS_FILE as a side effect.
    """
    # Ensure crossing list dependency is available/downloaded first
    _ = parse_crossing_list()

    if os.path.isfile(EPHEMERIS_FILE) and not force_rebuild:
        print("Ephemeris table already exists. Pass force_rebuild=True to regenerate.")
        return

    # Build a 60-second-cadence time grid spanning the full mission
    time_range = time_array(mission_start, mission_end, 60)

    ephemeris_data = parse_spice(time_range, units='km', frame=frame)

    # Save the table (overwrite only if force_rebuild, to avoid accidentally clobbering existing data otherwise)
    if force_rebuild:
        ephemeris_data.write(EPHEMERIS_FILE, overwrite = True)
    else:
        ephemeris_data.write(EPHEMERIS_FILE)

    print(f'Data Saved')
