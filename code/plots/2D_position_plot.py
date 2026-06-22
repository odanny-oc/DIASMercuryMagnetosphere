import numpy as np
from astropy.table import QTable
from sunpy.time import TimeRange

from hermpy.data import parse_messenger_mag
from hermpy.net import ClientMESSENGER
from hermpy.plotting import MultiPanel, SpectrogramPanel, TimeseriesPanel
import matplotlib.pyplot as plt


#Creates class to query MESSENGER data
c = ClientMESSENGER()

time_range = TimeRange("2011-09-27T12:00", "2011-09-27T17:00")

c.query(time_range, "MAG")
mag_file_paths = c.fetch()

mes_data: QTable = parse_messenger_mag(mag_file_paths, time_range)

pos_data = mes_data["X MSO", "Y MSO", "Z MSO"]
mag_data = mes_data["UTC", "Bx", "By", "Bz"]
time = mes_data["UTC"]

# Create arrays out of mag data to plot them later
Bx = np.array(mag_data["Bx"])
By = np.array(mag_data["By"])
Bz = np.array(mag_data["Bz"])

# Total magnetic field strength
mag_total = np.sqrt(mag_data["Bx"]**2 + mag_data["By"]**2 + mag_data["Bz"]**2)
mag_total_array = np.array(mag_total)

# Positional data
X = pos_data['X MSO']
Y = pos_data['Y MSO']
Z = pos_data['Z MSO']

# projected plot for cylindrical coordinates
P = np.sqrt(Y**2 +  Z**2)

# Mercury radius km
R = 2440

# Plot Mercury in 2D plane
t = np.linspace(0, 2*np.pi, 100)
X0 = R*np.cos(t)
X1 = R*np.sin(t)

fig, ax = plt.subplots(1,3)
#

def plane_plotter(X0, X1, Bx0, Bx1, plane, i):
    # Slicing of quiver data as to not make so noisy
    slicing = 100
    ax[i].plot(X0, X1, ls = '--', lw=1, color='r', label='MESSENGER trajectory')
    # Which plane are we plotting in
    ax[i].set_title(plane[0]+ '-' + plane[1] + ' Plane')
    # Plot magnetic field strength and direction in the X-Y plane
    ax[i].quiver(X0[::slicing], X1[::slicing], Bx0[::slicing], Bx1[::slicing], mag_total_array[::slicing], cmap='viridis', scale_units='width', scale=1000)
    ax[i].set_xlabel(plane[0] + r"$_{\text{MSO}}$ (km)")
    ax[i].set_ylabel(plane[1] + r"$_{\text{MSO}}$ (km)")

# Plot Mercury for each plot
for i in range(len(ax)):
    ax[i].plot(X0, X1, lw=2, color='k', label=r'Mercury') 
    ax[i].set_aspect('equal')

X = np.array(X)
Y = np.array(Y)
Z = np.array(Z)

plane_plotter(X, Y, Bx, By, ["X","Y"], 0)
plane_plotter(Y, Z, By, Bz, ["Y","Z"], 1)
plane_plotter(X, Z, Bx, Bz, ["X","Z"], 2)


# Line with colourmap based on magnetic field strength at that time and position
plt.colorbar(ax[2].quiver(X, Z, Bx, Bz, mag_total_array, cmap='viridis'), pad=0.1, label=r'Total Magnetic Field strength ($|B|$ nT)')


plt.legend()
plt.savefig("../data/2D_positional_data_mag.svg")
plt.show()
