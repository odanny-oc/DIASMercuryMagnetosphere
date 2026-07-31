import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import datetime as dt

from astropy.table import QTable
from astropy.time import Time 

from hermpymod.functions.ephemeris_downsampler import parse_apoapsis_data, parse_periapsis_data, parse_crossing_list
from hermpymod.functions.downsampled_positional_data import parse_spice_downsampled
from hermpymod.functions.encounters import parse_encounters_list
from hermpymod.functions.plot_all import plot_all_ephemeris
from hermpymod.classes.panels import HistogramPanel


# Load crossing and encounter data set for plots
crossing_data = parse_crossing_list()
encounters_data = parse_encounters_list()


# Load apoapsis and periapsis data
periapsis_data = parse_periapsis_data()
periapsis_times = Time(periapsis_data["UTC"]).to_datetime()

apoapsis_data = parse_apoapsis_data()
apoapsis_times = Time(apoapsis_data["UTC"]).to_datetime()

# Mask out crossings by type BS or MP
bs_crossings = QTable([i for i in crossing_data if "BS" in i["Label"]])
bs_crossings_times = Time(bs_crossings["UTC"]).to_datetime()

mp_crossings = QTable([i for i in crossing_data if "MP" in i["Label"]])
mp_crossings_times = Time(mp_crossings["UTC"]).to_datetime()


# Create lists to break up by orbit. We track the nearest apoapsis or periapsis and the proceeding apoapsis and periapsis
dt_bs_next_apo_12_hour = []
dt_bs_next_apo_8_hour = []
dt_mp_next_peri_12_hour = []
dt_mp_next_peri_8_hour = []

# Last 12 hour orbit 800
# First 8 hour orbit 812

dt_bs_apoapsis_12_hour = []
dt_bs_apoapsis_8_hour = []
dt_mp_periapsis_12_hour = []
dt_mp_periapsis_8_hour = []

mp_shortest = []
bs_longest = []


for bs in bs_crossings_times:

    # Find nearest apoapsis ahead of given bs crossing
    i0 = np.searchsorted(apoapsis_times, bs)
    dt0 = apoapsis_times[i0] - bs

    # Add that to the relevant proceeding apoapsis list
    if i0 <=800:
        if dt0.total_seconds()/3600 >= 12:
            bs_longest.append([i0, bs])
        dt_bs_next_apo_12_hour.append(dt0)
    elif i0 >=812:
        if dt0.total_seconds()/3600 >= 8:
            bs_longest.append([i0, bs])
        dt_bs_next_apo_8_hour.append(dt0)
    else:
        pass

    # Check to see if previous apoapsis is closer
    dt1 = apoapsis_times[i0-1] - bs

    dts = np.array([abs(dt0), abs(dt1)])
    min_dt = min(dts)
    indices = [i0, i0 - 1]

    index = indices[np.argmin(dts)]

    
    if index <= 800:
        dt_bs_apoapsis_12_hour.append(min_dt)
    elif index >= 812:
        dt_bs_apoapsis_8_hour.append(min_dt)
    else:
        continue


#Do the same for magenetopause crossings and the periapsides

for mp in mp_crossings_times:
    i0 = np.searchsorted(periapsis_times, mp)
    dt0 = periapsis_times[i0] - mp

    if i0 <= 800:
        dt_mp_next_peri_12_hour.append(abs(dt0))
    elif i0 >= 812:
        dt_mp_next_peri_8_hour.append(abs(dt0))
    else:
        pass

    dt1 = periapsis_times[i0-1] - mp

    dts = np.array([abs(dt0), abs(dt1)])
    min_dt = min(dts)
    indices = [i0, i0 - 1]

    index = indices[np.argmin(dts)]

    # Save times that seem too close to periapsis
    if min_dt.total_seconds()/3600 <= 0.167:
        mp_shortest.append([index, mp])
    
    if index <= 800:
        dt_mp_periapsis_12_hour.append(min_dt)
    elif index >= 812:
        dt_mp_periapsis_8_hour.append(min_dt)
    else:
        continue


# Plots

plot_configs_mp_bs = [
        (dt_bs_apoapsis_12_hour, dt_bs_apoapsis_8_hour, "Bow shock crossings to nearest apoapsis", "yellow", "orange", np.arange(0,6, 0.167)),
        (dt_mp_periapsis_12_hour, dt_mp_periapsis_8_hour ,"Magenetopause crossings to nearest periapsis", "purple", "pink", np.arange(0,2, 0.167)),

        (dt_bs_next_apo_12_hour, dt_bs_next_apo_8_hour, "Bow shock crossings to next apoapsis", "yellow", "orange", np.arange(0,13, 0.167)),
        (dt_mp_next_peri_12_hour, dt_mp_next_peri_8_hour, "Magenetopause crossings to next periapsis", "purple", "pink", np.arange(0,12, 0.167))
        ]

"""
Plot each data set together with both orbits and different colors, and seperate in two different subplots
"""

for data_12, data_8, label, color1, color2, bins in plot_configs_mp_bs:
    # Change data to numerics
    data_12 = [i.total_seconds()/3600 for i in data_12]
    data_8 = [i.total_seconds()/3600 for i in data_8]
    
    # Make histograms

    hist_12 = HistogramPanel(data_12, bins=bins, color=color1, minmax=True, average=False, label=f"12 Hour Orbits \nNumber of points {len(data_12)}",zorder=3)
    hist_8 = HistogramPanel(data_8, bins=bins, color=color2, minmax=True, average=False, label=f"8 Hour Orbits \nNumber of points {len(data_8)}", zorder=1)

    hists = [hist_12, hist_8]

    # Label and title histograms
    for hist in hists:
        hist.ax_set_params = {
                "title": f"{label}",
                "xlabel": "Time (hours)",
                "ylabel": "Number of crossings",
                "yscale": "log"
                }
    # Plot in different subplots
    hist = hist_12 + hist_8
    hist.plot(show=False)

    # Plot one on top of other (12 hour in front of 8 hour as it has less frequency i.e fewer 12 hour orbits)
    fig, ax = plt.subplots()
    for hist in hists:
        hist._plot_on(ax)

    # Labels
    plot_handles, _ = ax.get_legend_handles_labels()
    ax.set_xlabel("Time (Hours)")
    ax.set_ylabel("Number of crossings")
    ax.set_yscale("log")
    fig.suptitle(label)
    handles = [Line2D([], [], color = 'none', label = f"Total number of points = {len(data_12) + len(data_8)}")]
    ax.legend(handles=plot_handles + handles)


"""
Analysis of crossings that are too close or far away
"""

extrema_config = [
        (bs_longest, apoapsis_times, "apoapsis", "Longest Times from BS to apoapsis"),
        (mp_shortest, periapsis_times, "periapsis", "Shortest Times from MP to periapsis")
        ]

# Take four examples over 8 hour time window
tw = dt.timedelta(hours=4)

for extrema, location_times, location, title in extrema_config:
    print(title)
    print(len(extrema))
    slicing = len(extrema) // 3
    
    times_dt = []

    for orbit, crossing in extrema:
        time_dt = (location_times[orbit] - crossing).total_seconds()/3600
        times_dt.append(time_dt)
        
    max_index = np.where(times_dt == np.max(times_dt))[0][0]
    print(np.where(times_dt == np.max(times_dt)))
    print(np.max(times_dt))
    print(extrema[max_index][1].isoformat(), location_times[extrema[max_index][0]])

    for orbit, crossing in extrema[::slicing]:
        print("Crossing " + crossing.isoformat())
        print("Orbit number " + str(orbit))
        print("Time of next " + location + " " + location_times[orbit].isoformat())
        # Define time range of plot
        time_range = [location_times[orbit] - tw, location_times[orbit] + tw]
        # Find time between crossing and apoapsis/periapsis
        time_dt = (location_times[orbit] - crossing).total_seconds()/60

        # Find if crossing was before or after apo/periapsis
        if time_dt < 0:
            direction = "after"
            time_range_crossing = [location_times[orbit], crossing]
        else:
            direction = "to"
            time_range_crossing = [crossing, location_times[orbit]]

        time_dt = abs(time_dt)

        # Plot all planes X-Y... and cylindrical coords of time range
        fig, ax = plot_all_ephemeris(time_range, crossings=crossing_data, encounters=encounters_data, add_legend=False)

        # Plot length till apo/periapsis on top of above plot
        plot_all_ephemeris(time_range_crossing, color='green', mercury=False, scatter=False, ax=ax, label=f"{time_dt:.2f} minutes {direction} {location}")
        ax[-1].set_title(title)
        handles, _ = ax[0].get_legend_handles_labels()
        fig.legend(handles=handles)

plt.show()
