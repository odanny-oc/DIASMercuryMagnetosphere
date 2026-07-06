import pandas as pd
import os
from astropy.time import Time


home_dir = os.getenv('HOME')
data_dir = os.path.join(home_dir , '.ephemeris_data/')

"""
Need to download Philpott list from https://borealisdata.ca/dataset.xhtml?persistentId=doi:10.5683/SP2/1U6FEO
"""

philpott_crossings = pd.read_csv(data_dir + 'philpott_list.tab', delimiter=\t)

data_gap_masks ={
        "9" : philpott_crossings["Boundary number"] != 9.0,
        "10" : philpott_crossings["Boundary number"] != 10.0
        }

philpott_crossings = philpott_crossings[data_gap_masks["9"]]
philpott_crossings = philpott_crossings[data_gap_masks["10"]]

boundarynumber_dict = {
        1 : 'BS_IN outer',
        2 : 'BS_IN inner',
        3 : 'MP_IN inner',
        4 : 'MP_IN outer',
        5 : 'MP_OUT inner',
        6 : 'MP_OUT outer',
        7 : 'BS_OUT inner',
        8 : 'BS_OUT outer',
        }

label = []
for i in philpott_crossings["Boundary number"]:
    label.append(boundarynumber_dict[i])
    

philpott_crossings["UTC"] = [
    f"{int(y)}:{int(d):03d}:{int(h):02d}:{int(m):02d}:{s:06.3f}"
    for y, d, h, m, s in zip(philpott_crossings["Year"], philpott_crossings["Day of year"], philpott_crossings["Hour"], philpott_crossings["Minute"], philpott_crossings["Second"])
    ]

philpott_dt = pd.DataFrame()

philpott_dt["UTC"] = philpott_crossings["UTC"]
philpott_dt["X MSO"] = philpott_crossings["X_MSO (km)"]
philpott_dt["Y MSO"] = philpott_crossings["Y_MSO (km)"]
philpott_dt["Z MSO"] = philpott_crossings["Z_MSO (km)"]
philpott_dt["Label"] = label

philpott_dt.to_csv(data_dir + 'philpott_crossings_list_2020.csv')

"""
Order Philpott data like Sun's encounter list
"""

boundarynumber_masks ={
        "1" : {"mask": philpott_crossings["Boundary number"] == 1.0, "time": "Time Start" , "label": "BSI"},
        "2" : {"mask": philpott_crossings["Boundary number"] == 2.0, "time": "Time End" , "label": "BSI"},
        "3" : {"mask": philpott_crossings["Boundary number"] == 3.0, "time": "Time Start", "label": "MPI"},
        "4" : {"mask": philpott_crossings["Boundary number"] == 4.0, "time": "Time End", "label": "MPI"},
        "5" : {"mask": philpott_crossings["Boundary number"] == 5.0, "time": "Time Start", "label": "MPO"},
        "6" : {"mask": philpott_crossings["Boundary number"] == 6.0, "time": "Time End" , "label": "MPO"},
        "7" : {"mask": philpott_crossings["Boundary number"] == 7.0, "time": "Time Start" , "label": "BSO"},
        "8" : {"mask": philpott_crossings["Boundary number"] == 8.0, "time": "Time End" , "label": "BSO"},
        }

rows = []
pairs = zip(["1","3","5","7"], ["2","4","6","8"])

for start_idx, end_idx in pairs:
    start_times = philpott_crossings["UTC"][boundarynumber_masks[start_idx]["mask"]]
    end_times   = philpott_crossings["UTC"][boundarynumber_masks[end_idx]["mask"]]
    label       = boundarynumber_masks[start_idx]["label"]
    
    for s, e in zip(start_times, end_times):
        rows.append({"Time Start": s, "Time End": e, "Label": label})

philpott_dt_encounter = pd.DataFrame(rows)
philpott_dt_encounter = philpott_dt_encounter.sort_values("Time Start").reset_index(drop=True)
philpott_dt_encounter["Time Start"] = [Time(i) for i in philpott_dt_encounter["Time Start"]]
philpott_dt_encounter["Time End"] = [Time(i) for i in philpott_dt_encounter["Time End"]]

philpott_dt_encounter["Time Start"] = pd.to_datetime([t.iso for t in philpott_dt_encounter["Time Start"]])
philpott_dt_encounter["Time End"] = pd.to_datetime([t.iso for t in philpott_dt_encounter["Time End"]])

encounter_duration = (philpott_dt_encounter["Time End"] - philpott_dt_encounter["Time Start"])

encounter_duration = [i.total_seconds()/3600 for i in encounter_duration]

philpott_dt_encounter.insert(loc=3, column="Encounter Duration", value=encounter_duration)

philpott_dt_encounter.to_csv(data_dir + 'philpott_encounter_list_2020.csv', index=False)
