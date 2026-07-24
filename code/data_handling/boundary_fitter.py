import numpy as np
import matplotlib.pyplot as plt
import datetime as dt

from astropy.time import Time
import astropy.units as u
from astropy.table import QTable
from hermpy.utils.constants import Constants
from hermpymod.functions.encounters import encounter_finder, parse_crossing_list 
from hermpymod.classes.panels import PlanarplotPanel
from hermpymod.functions.ephemeris_downsampler import abs_r
from hermpymod.functions.boundary_models_mod import plot_magnetospheric_boundaries
from hermpymod.functions.plot_all import plot_all_ephemeris


Zd = Constants.DIPOLE_OFFSET.to("Mercury Radii")


def shue_model(cos_theta, Rss, alpha):
    return Rss*(2/(1 + cos_theta))**alpha


crossing_data = parse_crossing_list()

encounter_data = encounter_finder(crossing_data)

boundary_points = []


for encounter in encounter_data:
    time = encounter["UTC"][0]
    x_pos = np.average(encounter["X MSO"])
    y_pos = np.average(encounter["Y MSO"])
    z_pos = np.average(encounter["Z MSO"])
    label = encounter["Label"][0][:2] + encounter["Trajectory Direction"][0][0]
    boundary_points.append([time,x_pos,y_pos,z_pos, label])

boundary_points = np.array(boundary_points)

positions =[boundary_points[:,1].astype(float), boundary_points[:,2].astype(float), boundary_points[:,3].astype(float)] 

distance = abs_r(positions)

average_crossing_data = QTable({
    "UTC": boundary_points[:,0],
    "|R|": distance * u.Unit("Mercury Radii"),
    "X MSO": positions[0] * u.Unit("Mercury Radii"),
    "Y MSO": positions[1] * u.Unit("Mercury Radii"),
    "Z MSO": positions[2] * u.Unit("Mercury Radii"),
    "Label": boundary_points[:,4],
    })

mp_crossings = [crossing for crossing in average_crossing_data if "MP" in crossing["Label"]]

boundary_parameter = []

for idx, mp in enumerate(mp_crossings):
    try:
        next = mp_crossings[idx + 1]
    except IndexError:
        continue
    if mp["Label"] != next["Label"]:
        rho1 = np.sqrt(mp["Y MSO"]**2 + (mp["Z MSO"] - Zd)**2)
        rho2 = np.sqrt(next["Y MSO"]**2 + (next["Z MSO"] - Zd)**2)

        r1 = np.sqrt(mp["X MSO"]**2 + rho1**2)
        r2 = np.sqrt(next["X MSO"]**2 + rho2**2)

        cos_theta1 = mp["X MSO"]/r1
        cos_theta2 = next["X MSO"]/r2

        alpha = np.log(r1/r2)/np.log(shue_model(cos_theta1,Rss=1,alpha=1)/shue_model(cos_theta2,Rss=1,alpha=1))
        r0 = r1/shue_model(cos_theta1,Rss=1, alpha=alpha)

        boundary_parameter.append([float(r0.value), float(alpha.value)])

tw = dt.timedelta(hours=4)

print(len(boundary_parameter))

boundary_parameter = np.array(boundary_parameter)

sub_solar_magnetopause_params = boundary_parameter[:,0]
alpha_params = boundary_parameter[:,1]


print(sub_solar_magnetopause_params.max())
max_index = np.where(sub_solar_magnetopause_params == sub_solar_magnetopause_params.max())[0][0]
second_largest = sorted(sub_solar_magnetopause_params)[-10:-1]
second_smallest = sorted(sub_solar_magnetopause_params)[:10]
# second_max_index = np.where(sub_solar_magnetopause_params==second_largest)[0][0]

print(second_largest, second_smallest)


print(np.average(boundary_parameter[:,0]), np.average(boundary_parameter[:,1]))

slicing = len(boundary_parameter)//3

plot_range = [i for i in range(len(boundary_parameter))[::slicing]]

plot_range.append(max_index)
print(len(plot_range))
plot_range = [max_index]

for i in plot_range:
    time_range = [average_crossing_data["UTC"][i].to_datetime() - tw, average_crossing_data["UTC"][i].to_datetime() + tw]

    fig, ax = plot_all_ephemeris(time_range, crossings=average_crossing_data)

    config = [
            (ax[0], "xy", False),
            (ax[1], "xz", False),
            (ax[2], "yz", False),
            (ax[3], "xy", True),
            ]
    for axis, plane, cylindrical in config:
        plot_magnetospheric_boundaries(axis, plane=plane, frame="MSO", alpha=boundary_parameter[i][1], sub_solar_magnetopause=boundary_parameter[i][0], color='red', cylindrical=cylindrical)

plt.show()
