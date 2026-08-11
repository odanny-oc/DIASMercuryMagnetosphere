from astropy.table import Table
import pandas as pd
import os
import subprocess
from pathlib import Path

os.chdir(Path(__file__).resolve().parent)

from hermpymod.paths import DATA_DIR
data_dir = DATA_DIR

"""
This file creates the Sun dataset it a more usable format
"""


"""
These are some functions to parse the .txt files that the Sun data is presentes as
"""

def parse_bsi_line(line):
    # 6 tokens for start, 6 for end, then "BSI", then qualifier
    t = line.split()
    if len(t) < 14:
        return None
    
    # Recreate strings for pd formatting
    def to_isot(parts):
        y, mo, d, h, mi, s = parts
        s = max(0.0, float(s))  # clamp negatives
        return f"{y}-{mo}-{d}T{h}:{mi}:{float(s):06.3f}"
    
    t_start = to_isot(t[0:6])
    t_end   = to_isot(t[6:12])
    label = t[12] # BSO MSI etc.
    qualifier = t[13] # m or s

    dict = {
        "Time Start": t_start,
        "Time End":   t_end,
        "Label": label,
        "Type":    qualifier,
    }

    return dict


def sun_data_parser(data):
    """
    Function to create pandas dataframe of the Sun data
    """
    rows = []
    with open(data) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parsed = parse_bsi_line(line)
            if parsed:
                rows.append(parsed)

    df = pd.DataFrame(rows)
    df["Time Start"] = pd.to_datetime(df["Time Start"])
    df["Time End"]   = pd.to_datetime(df["Time End"])
    df["Encounter Duration"] = (df["Time End"] - df["Time Start"]).dt.total_seconds()/3600
    return df

home_dir = os.getenv('HOME')
os.makedirs(data_dir, exist_ok = True)

"""
Below downloads the data and combines it into a time ordered CSV
"""

urls = ['https://zenodo.org/records/18236915/files/Bow_Shock_In_Time_Duration_public_version_WeijieSun_20260113.txt?download=1',
        'https://zenodo.org/records/18236915/files/Bow_Shock_Out_Time_Duration_public_version_WeijieSun_20260113.txt?download=1',
        'https://zenodo.org/records/8298647/files/MagPause_In_Time_Duration__public_version_WeijieSun_20230829.txt?download=1',
        'https://zenodo.org/records/8298647/files/MagPause_Out_Time_Duration_public_version_WeijieSun_20230829.txt?download=1']

tmp_dir = os.path.join(data_dir, "tmp/")

for url in urls:
    os.makedirs(tmp_dir, exist_ok=True)
    # Create name from URL
    tempname = url.split('/')[-1].split('?')[0]
    save_path = os.path.join(tmp_dir, tempname)
    subprocess.run(["curl", "-Lo", save_path, url])
    

data_frames = []
# Create dataframes
for file in os.listdir(tmp_dir):
    data_frames.append(sun_data_parser(tmp_dir + file))

sun_crossing_data = pd.concat(data_frames, ignore_index=True)
# Time order data
sun_crossing_data.sort_values("Time Start", inplace=True)

# Save data to CSV
sun_crossing_data.to_csv(os.path.join(data_dir, "sun_2023_crossing.csv"), index=False)
