import xarray as xr
import numpy as np
from astropy.table import QTable
from sunpy.time import TimeRange

from hermpy.data import parse_messenger_fips, parse_messenger_mag
from hermpy.net import ClientMESSENGER
from hermpy.plotting import MultiPanel, SpectrogramPanel, TimeseriesPanel
import matplotlib.pyplot as plt
from matplotlib import colormaps as cm
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib.colors import Normalize
from PIL import Image

# Load the texture
img = Image.open("../plots_and_images/mercury_texture.jpg")
img = np.array(img) / 255.0  # normalise to [0, 1]

#Creates class to query MESSENGER data
c = ClientMESSENGER()
time_range = TimeRange("2011-09-27T12:00", "2011-09-27T23:00")

c.query(time_range, "MAG")
mag_file_paths = c.fetch()

# %%
# We use Astropy's Quantity Table (QTable) for time series data, and Xarray's
mes_data: QTable = parse_messenger_mag(mag_file_paths, time_range)

#keep only spatial data, discarding the magnetic field strength
pos_data = mes_data["X MSO", "Y MSO", "Z MSO"]
pos_and_time_data = mes_data["UTC", "X MSO", "Y MSO", "Z MSO"]
mag_data = mes_data["UTC", "Bx", "By", "Bz"]

mag_total = np.sqrt(mag_data["Bx"]**2 + mag_data["By"]**2 + mag_data["Bz"]**2)

slicing = len(pos_and_time_data)//8

pos_and_time_data = pos_and_time_data[::slicing]

mag_panel = TimeseriesPanel(mag_data)
fig, ax = mag_panel.plot(show=False)
ax.plot(mag_data["UTC"].to_datetime(leap_second_strict="warn") , mag_total.value, color='k', label="|B|")
ax.legend()

X = pos_data['X MSO']
Y = pos_data['Y MSO']
Z = pos_data['Z MSO']

# projected plot for cylindrical coordinates
P = np.sqrt(Y**2 +  Z**2)

u = np.linspace(0, 2 * np.pi, img.shape[0])  # columns = longitude
v = np.linspace(0, np.pi, img.shape[1])       # rows = latitude

#km
R = 2440

X0 = R * np.outer(np.cos(u), np.sin(v))
X1 = R * np.outer(np.sin(u), np.sin(v))
X2 = R * np.outer(np.ones(np.size(u)), np.cos(v))

ax = plt.figure(figsize=(10,10)).add_subplot(projection='3d')
#

ax.plot_surface(X0, X1, X2, facecolors = img, rstride=1, cstride=1, zorder=0, label=r'Mercury')

for i in pos_and_time_data:
    ax.plot(i[1].value, i[2].value, i[3].value, color='red', zorder=5, markersize=0.5)
    ax.text(i[1].value, i[2].value, i[3].value, i[0], fontsize=10)

# Line with colourmap based on magnetic field strength at that time and position
# ax.plot(X, Y, Z)

X = np.array(X)
Y = np.array(Y)
Z = np.array(Z)

points = np.vstack([X, Y, Z]).T
segments = np.concatenate([points[:-1, None, :], points[1:, None, :]], axis=1)

mag_total_array = np.array(mag_total)

lc = Line3DCollection(segments, array=mag_total_array, cmap='viridis', linewidth=2)

ax.set_xlim(min(X.min() - R, R - X.min()), X.max() + R)
ax.set_ylim(min(Y.min() - R, R - Y.min()), Y.max() + R)
ax.set_zlim(min(Z.min() - R, R - Z.min()), Z.max() + R)

ax.set_xlabel(r"$X_{\text{MSO}}$ (km)")
ax.set_ylabel(r"$Y_{\text{MSO}}$ (km)")
ax.set_zlabel(r"$Z_{\text{MSO}}$ (km)")

ax.add_collection(lc)

plt.colorbar(lc, ax=ax, pad=0.1, label=r'Total Magnetic Field strength ($|B|$ nT)')

ax.set_aspect('equal')

plt.savefig("../data/positional_data_mag.svg")
plt.show()
