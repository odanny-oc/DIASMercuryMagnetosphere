import numpy as np
import matplotlib.pyplot as plt
from astropy.table import QTable
from sunpy.time import TimeRange
from astropy.time import Time
import pickle

from hermpy.data import parse_messenger_fips, parse_messenger_mag
from hermpy.net import ClientMESSENGER
from hermpy.plotting import MultiPanel, SpectrogramPanel, TimeseriesPanel

import matplotlib.pyplot as plt
from matplotlib import colormaps as cm
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from PIL import Image
import datetime as dt
from datetime import timedelta

import matplotlib

matplotlib.use('QtAgg')

with open('../data/orbits_data.pkl', 'rb') as f:
    orbits = pickle.load(f)

# Number of crossings per 12hrs 
crossing_number = [len(orbits["orbit_times"][i]) for i in range(len(orbits))]

print(crossing_number)

average = np.mean(crossing_number)

fig, ax = plt.subplots()

ax.tick_params(labelsize=12)
# Plot histogram of number of crosses per orbit
ax.hist(crossing_number, bins = 20, edgecolor= 'k')

ax.axvline(average, ls='--', color='r', label=f'Average number of crossings = {average:.2f}')

ax.set_title("Number of magnetic field boundary crossings per 12hr orbit, (2011-02-23 - 2012-03-31)", fontsize=24)
ax.set_xlabel("Number of crossings", fontsize=16)
ax.set_ylabel("Number of orbits", fontsize=16)
print('Plotted Histogram')
plot_handles, plot_labels = ax.get_legend_handles_labels()
handles = [Line2D([], [], color = 'none', label = f"Number of orbits = {len(orbits)}")]
ax.legend(handles= plot_handles + handles, fontsize=14)
plt.savefig('./crossings_histogram.svg')
plt.show()
