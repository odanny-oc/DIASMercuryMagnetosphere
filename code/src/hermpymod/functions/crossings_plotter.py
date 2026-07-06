import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from astropy.table import QTable, vstack, Column
import astropy.units as u
from astropy.time import Time
from astropy.table import QTable, hstack, vstack
from sunpy.time import TimeRange

from hermpy.plotting import Panel, TimeseriesPanel
from hermpy.data import parse_messenger_fips, parse_messenger_mag
from hermpy.net import ClientSPICE, ClientMESSENGER
from hermpy.utils import Constants as c

import spiceypy as spice

import datetime as dt

import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from hermpymod.functions.ephemeris_downsampler import parse_crossing_list


def plot_crossings(t_start, t_end, ax):
    crossing_list = parse_crossing_list()
    crossing_time = crossing_list["UTC"]
    mask = (crossing_time >= t_start) & (crossing_time <= t_end)
    crossing_list = crossing_list[mask]

    mask_dict = {
                "BS_OUT":{"mask": crossing_list["Label"] == "BS_OUT", "color":'yellow'},
                "BS_IN" :{ "mask": crossing_list["Label"] == "BS_IN", "color": 'red'},
                "MP_OUT":{ "mask": crossing_list["Label"] == "MP_OUT","color": 'purple'},
                "MP_IN" :{ "mask": crossing_list["Label"] == "MP_IN", "color": 'blue'},
                    }

    if not isinstance(ax, np.ndarray):
        ax = [ax]


    for ax in ax:
        for mask in mask_dict:
            ax.vlines(crossing_list["UTC"][mask_dict[mask]["mask"]].to_datetime(), ymin=0, ymax=1, transform=ax.get_xaxis_transform(), color= mask_dict[mask]["color"], label = mask, ls='--')

