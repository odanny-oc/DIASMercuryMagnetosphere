import matplotlib.pyplot as plt
import numpy as np

import astropy.units as u
from astropy.time import Time
from astropy.table import QTable, hstack, vstack
from sunpy.time import TimeRange

from hermpy.plotting import Panel, TimeseriesPanel
from hermpy.utils import Constants as c

import datetime as dt

from hermpymod.functions.ephemeris_downsampler import parse_crossing_list
from hermpymod.functions.encounters import parse_encounters_list


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
                "UNPHYSICAL (SW - MSp)" :{ "mask": crossing_list["Label"] == "UNPHYSICAL (SW -> MSp)" , "color": 'gray'},
                "UNPHYSICAL (MSp - SW)" :{ "mask": crossing_list["Label"] == "UNPHYSICAL (MSp -> SW)" , "color": 'brown'},
                    }

    if not isinstance(ax, np.ndarray):
        ax = [ax]


    for i in range(len(ax)):
        for mask in mask_dict:
            if len(ax) == 1:
                ax[i].vlines(crossing_list["UTC"][mask_dict[mask]["mask"]].to_datetime(), ymin=0, ymax=1, transform=ax[i].get_xaxis_transform(), color= mask_dict[mask]["color"], label = mask, ls='--')
            else:
                ax[i].vlines(crossing_list["UTC"][mask_dict[mask]["mask"]].to_datetime(), ymin=0, ymax=1, transform=ax[i].get_xaxis_transform(), color= mask_dict[mask]["color"], label = mask if i == 1 else None, ls='--')


def plot_encounters(t_start, t_end, ax):
    encounter_list = parse_encounters_list()
    start_times = Time(encounter_list["Time Start"]).to_datetime()
    end_times = Time(encounter_list["Time End"]).to_datetime()
    start_mask = (start_times >= t_start) & (start_times <= t_end)
    end_mask = (end_times >= t_start) & (end_times <= t_end)
    mask = start_mask | end_mask
    encounter_list = encounter_list[mask]

    if not isinstance(ax, np.ndarray):
        ax = [ax]


    for i in range(len(ax)):
        if len(ax) == 1:
            for idx, encounter in enumerate(encounter_list):
                encounter_times = [Time(encounter["Time Start"]).to_datetime(), Time(encounter["Time End"]).to_datetime()]
                ax[i].axvspan(encounter_times[0], encounter_times[-1], alpha=0.7, color='orange', label=f"Encounter \nNumber of encounters {len(encounter_list)}")
        else:
            for idx, encounter in enumerate(encounter_list):
                encounter_times = [Time(encounter["Time Start"]).to_datetime(), Time(encounter["Time End"]).to_datetime()]
                ax[i].axvspan(encounter_times[0], encounter_times[-1], alpha=0.7, color='orange', label=f"Encounter \nNumber of encounters {len(encounter_list)}" if idx ==0 and i==1 else None)
