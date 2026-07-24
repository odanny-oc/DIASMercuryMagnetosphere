import matplotlib.pyplot as plt
import numpy as np

from astropy.table import QTable, vstack

from hermpymod.classes.panels import parse_crossing_list, plot_magnetospheric_boundaries
from hermpymod.functions.encounters import encounter_finder


crossing_data = parse_crossing_list()

mp_crossings = [crossing for crossing in crossing_data if "MP" in crossing["Label"]]
mp_crossings = vstack(mp_crossings)

mp_encounters = encounter_finder(mp_crossings)

mp_x_errors = []
mp_rho_errors = []


for encounter in mp_encounters:
    x_err = [abs(encounter[0]["X MSO"] - encounter[-1]["X MSO"])]
    rho_err = [abs(np.sqrt(encounter[0]["Y MSO"]**2 +  encounter[0]["Z MSO"]**2) - np.sqrt(encounter[-1]["Y MSO"]**2 +  encounter[-1]["Z MSO"]**2))]
    
    x_err.extend([0] * (len(encounter)-1))
    rho_err.extend([0] * (len(encounter)-1))
    mp_x_errors.append(x_err)
    mp_rho_errors.append(rho_err)


mp_x_errors = np.concatenate(mp_x_errors)
mp_rho_errors = np.concatenate(mp_rho_errors)

mp_rho = [np.sqrt(i["Y MSO"]**2 + i["Z MSO"]**2) for i in mp_crossings]

mp_phi = [np.arctan(abs(i["Z MSO"]/i["Y MSO"]))* 180/np.pi for i in mp_crossings]

fig, ax = plt.subplots()

sc = ax.scatter(mp_crossings["X MSO"], mp_rho, c=mp_phi, cmap='viridis', s=0.1)
colors = sc.to_rgba(mp_phi)


for x, rho, xerr, yerr, color in zip(mp_crossings["X MSO"], mp_rho, mp_x_errors, mp_rho_errors, colors):
    ax.errorbar(x, rho, fmt='none', xerr=xerr, yerr=yerr, color=color, zorder=0, marker=None)


plot_magnetospheric_boundaries(ax, frame='MSO', sub_solar_magnetopause=1.4, cylindrical=True)

fig.colorbar(sc, ax=ax, label=r'Azimuthal angle ($^\circ$)')

ax.set_title("All MP Crossings for MESSENGER Mission")
ax.set_xlabel("X MSO Mercury Radii")
ax.set_ylabel(r"$\rho$ MSO Mercury Radii")

ax.set_xlim(-5,3)
ax.set_ylim(0,4)

ax.legend()

plt.show()
