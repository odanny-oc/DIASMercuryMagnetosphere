import matplotlib.pyplot as plt
import numpy as np

from astropy.table import QTable, vstack

from hermpymod.classes.panels import parse_crossing_list, plot_magnetospheric_boundaries
from hermpymod.functions.encounters import encounter_finder


crossing_data = parse_crossing_list()

bs_crossings = [crossing for crossing in crossing_data if "BS" in crossing["Label"]]
bs_crossings = vstack(bs_crossings)

bs_encounters = encounter_finder(bs_crossings)

bs_x_errors = []
bs_rho_errors = []


for encounter in bs_encounters:
    x_err = [abs(encounter[0]["X MSO"] - encounter[-1]["X MSO"])]
    rho_err = [abs(np.sqrt(encounter[0]["Y MSO"]**2 +  encounter[0]["Z MSO"]**2) - np.sqrt(encounter[-1]["Y MSO"]**2 +  encounter[-1]["Z MSO"]**2))]
    
    x_err.extend([0] * (len(encounter)-1))
    rho_err.extend([0] * (len(encounter)-1))
    bs_x_errors.append(x_err)
    bs_rho_errors.append(rho_err)


bs_x_errors = np.concatenate(bs_x_errors)
bs_rho_errors = np.concatenate(bs_rho_errors)

bs_rho = [np.sqrt(i["Y MSO"]**2 + i["Z MSO"]**2) for i in bs_crossings]

bs_phi = [np.arctan(abs(i["Z MSO"]/i["Y MSO"]))* 180/np.pi for i in bs_crossings]

fig, ax = plt.subplots()

sc = ax.scatter(bs_crossings["X MSO"], bs_rho, c=bs_phi, cmap='viridis', s=0.1)
colors = sc.to_rgba(bs_phi)


for x, rho, xerr, yerr, color in zip(bs_crossings["X MSO"], bs_rho, bs_x_errors, bs_rho_errors, colors):
    ax.errorbar(x, rho, fmt='none', xerr=xerr, yerr=yerr, color=color, zorder=0, marker=None)


plot_magnetospheric_boundaries(ax, frame='MSO', sub_solar_magnetopause=1.4, cylindrical=True, add_legend=True)

fig.colorbar(sc, ax=ax, label=r'Azimuthal angle ($^\circ$)')

ax.set_title("All BS Crossings for MESSENGER Mission")
ax.set_xlabel("X MSO Mercury Radii")
ax.set_ylabel(r"$\rho$ MSO Mercury Radii")
ax.set_xlim(-8,6)
ax.set_ylim(0,10)

ax.legend()

plt.show()
