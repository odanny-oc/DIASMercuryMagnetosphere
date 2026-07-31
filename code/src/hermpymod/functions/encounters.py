from astropy.table import QTable, vstack
from astropy.time import Time

import datetime as dt
import os
from hermpymod.functions.ephemeris_downsampler import parse_crossing_list


home_dir = os.getenv("HOME")
data_dir = os.path.join(home_dir, ".ephemeris_data/")
os.makedirs(data_dir, exist_ok = True)

crossing_data = parse_crossing_list()

"""
From list of crossings with direction, construct encounters list
"""
def encounter_finder(crossing_data=crossing_data):
    crossing_time = Time(crossing_data["UTC"]).to_datetime()
    crossing_label = crossing_data["Label"]
    crossing_direction = crossing_data["Trajectory Direction"]
    encounters_data = []
    cross_number = []

    for i in range(len(crossing_time)):
        crossing_type = crossing_label[i][:2]

        try:
            time_delta = crossing_time[i+1] - crossing_time[i]
            if crossing_type in crossing_label[i+1] and crossing_direction[i] == crossing_direction[i+1] and time_delta<=dt.timedelta(hours=13):
                cross_number.append(crossing_data[i])

            else:
                cross_number.append(crossing_data[i])
                rows = vstack(cross_number)

                encounters_data.append(rows)
                cross_number = []

        except IndexError:
            cross_number.append(crossing_data[i])
            rows = vstack(cross_number)

            encounters_data.append(rows)
            cross_number = []

    return encounters_data


"""
Creates/loads list of encounter times and durations
"""

def parse_encounters_list(force_rebuild=False):

    encounters_list_dir = os.path.join(data_dir, 'hollman_encounters_list_2025.csv')


    if force_rebuild:
        try:
            os.remove(encounters_list_dir)
        except FileNotFoundError:
            pass

    # Get crossings if not already downloaded
    try:
        encounters_data = QTable.read(encounters_list_dir)
    except FileNotFoundError:
        print("File not found. Building encounter list")

        encounters_list = encounter_finder()

        time_start = [i["UTC"][0] for i in encounters_list]
        time_end = [i["UTC"][-1] for i in encounters_list]

        label = []

        for encounter in encounters_list:
            encounter_type = encounter["Label"][0][:2]
            direction = encounter["Trajectory Direction"][0][0].upper()
            label.append(encounter_type + direction)

        encounter_duration = [(Time(time_end[i]).to_datetime() - Time(time_start[i]).to_datetime()).total_seconds()/3600 for i in range(len(time_start))]

        encounters_data = QTable({
                "Time Start": time_start,
                "Time End": time_end,
                "Label": label,
                "Encounter Duration": encounter_duration
                })
        encounters_data.write(encounters_list_dir)

    return encounters_data
