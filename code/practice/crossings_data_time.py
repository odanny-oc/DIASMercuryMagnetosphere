import numpy as np
import matplotlib.pyplot as plt
from astropy.table import QTable
from numpy._core import long
from sunpy.time import TimeRange
from astropy.time import Time
import pickle

from hermpy.data import parse_messenger_fips, parse_messenger_mag
from hermpy.net import ClientMESSENGER
from hermpy.plotting import MultiPanel, SpectrogramPanel, TimeseriesPanel

import matplotlib.pyplot as plt
from matplotlib import colormaps as cm
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib.colors import Normalize
from PIL import Image
import datetime as dt
from datetime import timedelta

# Load list of crossing ref Hollman 2025
crossing_list = QTable.read("../data/hollman_2025_crossing_list.csv")

# Convert from string to datetime
crossing_list["Time"] = Time(crossing_list["Time"]).to_datetime()

time = crossing_list["Time"]

# Cutoff before MESSENGER goes to 8 hr orbit
cutoff= dt.datetime(2012, 4, 1)

long_orbit = [i for i in time if i <= cutoff]

orbits = {}
crossing_type = {}

# MESSENGER orbit duration
dt_orbit = timedelta(hours=12)

starttime = long_orbit[0]

# floor divide, only includes full orbits
number_of_orbits = ((long_orbit[-1] - long_orbit[0])//dt_orbit)

orbit_number1 = 0
orbit_number2 = 1


# Divide each 12hr period into its own list
for i in range(number_of_orbits):
    print(i)
    # Take each point that lies between ith 12th hour period
    orbits[i] = [dates for dates in long_orbit if dates > starttime + orbit_number1 * dt_orbit and dates < starttime + orbit_number2 * dt_orbit]

    # Consider empty case
    if orbits[i] == []:
        crossing_type[i] = []
    else:
        # Fill crossing type list with matching index
        crossing_type[i] = crossing_list["Label"][long_orbit.index(orbits[i][0]):long_orbit.index(orbits[i][-1])+1]
    # Iterate the orbit
    orbit_number1 += 1
    orbit_number2 += 1

orbits_with_crossing_type = QTable({
    "orbit_times": list(orbits.values()),
    "crossing_types": list(crossing_type.values())
    })


# Save the data
with open('../data/orbits_data.pkl', 'wb') as f:
    pickle.dump(orbits_with_crossing_type, f)
