from astropy.table import QTable, vstack
import datetime as dt

def encounter_finder(crossing_data):
    crossing_time = crossing_data["UTC"].to_datetime()
    crossing_label = crossing_data["Label"]
    crossing_direction = crossing_data["Trajectory Direction"]
    encounters_data = []
    cross_number = []

    for i in range(len(crossing_time) - 1):
        time_delta = crossing_time[i+1] - crossing_time[i]
        crossing_type = crossing_label[i][:2]

        if crossing_type in crossing_label[i+1] and crossing_direction[i] == crossing_direction[i+1] and time_delta<=dt.timedelta(hours=13):
            cross_number.append(crossing_data[i])

        else:
            cross_number.append(crossing_data[i])
            rows = vstack(cross_number)

            encounters_data.append(rows)
            cross_number = []

    return encounters_data
