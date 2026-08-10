import os
from pathlib import Path
os.chdir(Path(__file__).resolve().parent)

import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import ks_2samp, wasserstein_distance

from hermpy.plotting import Panel, TimeseriesPanel
from hermpy.utils import Constants as c
import spiceypy as spice

from astropy.table import QTable, hstack, vstack, Column
import astropy.units as u
from astropy.time import Time

import os

from hermpymod.classes.panels import PlanarplotPanel, HistogramPanel
from hermpymod.functions.encounters import parse_encounters_list


home_dir = os.getenv('HOME')
from hermpymod.paths import DATA_DIR
data_dir = DATA_DIR
img_dir = "../../../plots_and_images/"

plt.style.use(img_dir + "presentation.mplstyle")


def dt_encounters(data, mode='same'):
    start_times = Time(data["Time Start"]).to_datetime()
    end_times = Time(data["Time End"]).to_datetime()
    dt_mp = []
    dt_bs = []
    check = (lambda a, b: a in b) if mode == 'same' else (lambda a, b: a not in b)

    for i in range(len(data) - 1):
        label = data[i]['Label'][:2]
        if check(label, data['Label'][i + 1]):
            delta_time = start_times[i+1] - end_times[i]
            # if delta_time > dt.timedelta(hours=13):
            #     continue
            if label == 'MP':
                dt_mp.append(delta_time)
            else:
                dt_bs.append(delta_time)
        else:
            continue
    dt_mp = [i.total_seconds()/3600 for i in dt_mp]
    dt_bs = [i.total_seconds()/3600 for i in dt_bs]
    return dt_bs, dt_mp


def cumdis(data):
    X = np.sort(data)
    F = np.array(range(len(data)))/float(len(data))
    return X,F


class CDFPanel(Panel):
    def __init__(self, data, label):
        super().__init__()
        if isinstance(data, list):
            self._data = data
            self._label = label
            self.k2D, self.k2p = ks_2samp(data[0], data[1])
            self.wasserdist = wasserstein_distance(data[0], data[1])
        else:
            self._data = [data]
            self._label = [label]
        self.cdf = [cumdis(i) for i in self._data]

    def _plot_on(self, ax):
        for i in range(len(self._data)):
            ax.plot(self.cdf[i][0], self.cdf[i][1], label=self._label[i])
        handles, _ = ax.get_legend_handles_labels()
        if len(self._data) == 2:
            legend_elements = [
                    Line2D([], [], color='none', label = f"D = {self.k2D:.2f}\np_value = {self.k2p:.2f}\nWasserstein Distance = {self.wasserdist:.2f}")
                    ]
            handles += legend_elements
        ax.set_ylabel("CDF")
        ax.grid()
        ax.legend(handles=handles)


phil_encounter_data = QTable.read(data_dir + 'philpott_encounter_list_2020.csv')

sun_encounter_data = QTable.read(data_dir + 'sun_2023_crossing.csv')

hollman_encounter_data = parse_encounters_list()


phil_dt_bs, phil_dt_mp = dt_encounters(phil_encounter_data)

sun_dt_bs, sun_dt_mp = dt_encounters(sun_encounter_data)

holl_dt_bs, holl_dt_mp = dt_encounters(hollman_encounter_data)

bins = np.arange(0, max(sun_dt_bs) + 2, 0.167)

configs = [
        (hollman_encounter_data, "Hollman"),
        (sun_encounter_data, "Sun et al."),
        (phil_encounter_data, "Philpott et al."),
        ]
        
for data, name in configs:

    dts = []
    for i in range(len(data) - 1):
        dts.append(Time(data["Time Start"][i+1]).to_datetime() - Time(data["Time End"][i]).to_datetime())


    dts = [i.total_seconds()/3600 for i in dts]
    bs, mp = dt_encounters(data)
    bs_mp = [
            (bs, "BS", "yellow"),
            (mp, "MP", "purple"),
            (dts, "all", "Green"),
            ]
    hists = []
    for data_type, label, color in bs_mp:
        hist = HistogramPanel(data_type, bins = bins, color=color)
        hist.ax_set_params = {
                "title": f"Time between {label} Encounters ({name})",
                "xlabel": "Time (Hours)",
                "ylabel": "Number of crossings",
                "yscale": "log",
                "xlim" : (0,12)
                }
        hists.append(hist)

    plot = hists[0] + hists[1]
    plot.plot(show=False, figsize=(18,16))
    plt.savefig(img_dir + f"{name}_dt_encounters_bs_mp.svg")
    hists[-1].plot(show=False)
    plt.savefig(img_dir + f"{name}_dt_encounters_all.svg")


    encounter_hist = HistogramPanel(data["Encounter Duration"], bins=bins)
    encounter_hist.ax_set_params = {
            "title" : f"Duration of Encounters ({name})",
            "xlabel": "Time (hours)",
            "ylabel": "Number of Encounters",
            "yscale": "log",
            "xlim" : (0,8)
            } 
    encounter_hist.plot(show=False)
    plt.savefig(img_dir + f"{name}_dt_encounters.svg")


bs_cdf_sun_phil = CDFPanel([sun_dt_bs, phil_dt_bs], label=["Sun", "Philpott"])
bs_cdf_holl_phil = CDFPanel([holl_dt_bs, phil_dt_bs], label=["Hollman", "Philpott"])
bs_cdf_holl_sun = CDFPanel([holl_dt_bs, sun_dt_bs], label=["Hollman", "Sun"])

bs_cdf_list = [bs_cdf_sun_phil, bs_cdf_holl_phil, bs_cdf_holl_sun]


for i in bs_cdf_list:
    i.ax_set_params = {
            "title" :"Cumulative Distribtution Functions (CDF) for Time between BS Encounters",
            "xlabel": "Time (Hours)",
            "xlim": (0,15),
    }


mp_cdf_sun_phil = CDFPanel([sun_dt_mp, phil_dt_mp], label=["Sun", "Philpott"])
mp_cdf_holl_phil = CDFPanel([holl_dt_mp, phil_dt_mp], label=["Hollman", "Philpott"])
mp_cdf_holl_sun = CDFPanel([holl_dt_mp, sun_dt_mp], label=["Hollman", "Sun"])

mp_cdf_list = [mp_cdf_sun_phil, mp_cdf_holl_phil, mp_cdf_holl_sun]

for i in mp_cdf_list:
    i.ax_set_params = {
            "title" :"Cumulative Distribtution Functions (CDF) for Time between MP Encounters",
            "xlabel": "Time (Hours)",
            "xlim": (0,11),
    }

cdf_sun_phil = bs_cdf_sun_phil + mp_cdf_sun_phil
cdf_holl_phil = bs_cdf_holl_phil + mp_cdf_holl_phil
cdf_holl_sun = bs_cdf_holl_sun + mp_cdf_holl_sun

cdf_sun_phil.plot(show=False, sharex=False, figsize=(18,16))
plt.savefig(img_dir + "cdf_sun_phil.svg")
cdf_holl_phil.plot(show=False, sharex=False, figsize=(18,16))
plt.savefig(img_dir + "cdf_hollman_phil.svg")
cdf_holl_sun.plot(show=False, sharex=False, figsize=(18,16))
plt.savefig(img_dir + "cdf_hollman_sun.svg")

for data, name in configs:
    bs, mp = dt_encounters(data, mode='diff')
    bs_mp = [
            (bs, "BS-MP", "red"),
            (mp, "MP-BS", "blue"),
            ]
    hists = []
    for data_type, label, color in bs_mp:
        hist = HistogramPanel(data_type, bins = bins, color=color)
        hist.ax_set_params = {
                "title": f"Time between {label} Encounters ({name})",
                "xlabel": "Time (Hours)",
                "ylabel": "Number of encounters",
                "yscale": "log",
                "xlim" : (0,12)
                }
        hists.append(hist)

    plot = hists[0] + hists[1]
    plot.plot(show=False, figsize=(18,16))


plt.show()
