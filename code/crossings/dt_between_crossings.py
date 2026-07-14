import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import QTable, vstack, Column
from hermpy.plotting import Panel, TimeseriesPanel
from hermpy.net import ClientSPICE, ClientMESSENGER
from astropy.time import Time
import spiceypy as spice
from hermpy.utils import Constants as c
import datetime as dt
from matplotlib.lines import Line2D
from astropy.table import QTable, hstack, vstack
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from hermpymod.classes.panels import PlanarplotPanel
from hermpymod.classes.panels import HistogramPanel
from hermpymod.functions.ephemeris_downsampler import parse_crossing_list


home_dir = os.getenv('HOME')
data_dir = os.path.join(home_dir, '.ephemeris_data/')
img_dir = "../../plots_and_images/"

plt.style.use(img_dir + "presentation.mplstyle")

hollman_crossing_list = parse_crossing_list()


def delta_crossing_times(crossing_list):
    crossing_time = Time(crossing_list["UTC"]).to_datetime()
    crossing_label = crossing_list["Label"]

    crossings_delta_t = [(crossing_time[i+1] - crossing_time[i]) for i in range(len(crossing_time) -1)]

    mp_mask = ["MP" in i for i in crossing_label]
    bs_mask = ["BS" in i for i in crossing_label]

    crossing_time_mp = crossing_time[mp_mask]
    crossing_time_bs = crossing_time[bs_mask]

    crossings_delta_t = [i.total_seconds()/3600 for i in crossings_delta_t]

    crossings_delta_t_mp = [(crossing_time_mp[i+1] - crossing_time_mp[i]) for i in range(len(crossing_time_mp) -1)]

    crossings_delta_t_bs = [(crossing_time_bs[i+1] - crossing_time_bs[i]) for i in range(len(crossing_time_bs) -1)]

    crossings_delta_t_mp = [i.total_seconds()/3600 for i in crossings_delta_t_mp]
    crossings_delta_t_bs = [i.total_seconds()/3600 for i in crossings_delta_t_bs]
    return crossings_delta_t, crossings_delta_t_bs, crossings_delta_t_mp

configs = [
        (hollman_crossing_list, "Hollman")
        ]

bins = np.arange(0,10,0.167)
for data, name in configs:
    dt, bs, mp = delta_crossing_times(data)

    data_sets=[
            (dt, "all", "C0"),
            (bs, "BS", "yellow"),
            (mp, "MP", "purple"),
            ]

    plots = []

    for data_type, label, color in data_sets:
        hist = HistogramPanel(data_type, bins = bins, color=color)
        hist.ax_set_params = {
                "title": f"Time between {label} crossings",
                "xlabel": "Time (Hours)",
                "ylabel": "Number of crossings",
                "yscale": "log"
                }
        plots.append(hist)

    dt_hist_bs_mp = plots[1] + plots[2]
    dt_hist = plots[0]
    dt_hist.plot(show=False)
    plt.savefig(img_dir + "hollman_crossings_dt.svg")
    dt_hist_bs_mp.plot(show=False, figsize=(18,16))
    plt.savefig(img_dir + "hollman_crossings_dt_bs_mp.svg")

plt.show()
