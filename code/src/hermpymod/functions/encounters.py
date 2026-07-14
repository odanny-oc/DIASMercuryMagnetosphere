from astropy.table import QTable, vstack
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
    crossing_time = crossing_data["UTC"].to_datetime()
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

    save_path = os.path.join(data_dir, "hollman_encounters_list_2025.csv")


    if force_rebuild:
        os.remove(save_path)
        os.remove(encounters_list_dir)

    # Get crossings if not already downloaded
    try:
        encounters_data = QTable.read(encounters_list_dir)
    except FileNotFoundError:
        crossing_time = Time(crossing_data["UTC"]).to_datetime()
        crossing_label = crossing_data["Label"]
        crossing_direction = crossing_data["Trajectory Direction"]
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
                cross_number.append(crossing_data[i])
                dts.append(time_delta)
            else:
                cross_number.append(crossing_data[i])
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

        encounters_data = QTable({
                "Time Start": time_start,
                "Time End": time_end,
                "Label": label,
                "Encounter Duration": encounter_duration
                })
        encounters_data.write(data_dir + "hollman_encounters_list_2025.csv")

    return encounters_data
