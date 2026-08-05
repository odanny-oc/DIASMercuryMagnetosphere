import numpy as np
from astropy.table import QTable
from astropy.time import Time

import astropy.units as u

import pickle
import spiceypy as spice

from hermpy.net import ClientSPICE

import datetime as dt
from datetime import timedelta

import os

from hermpymod.functions.ephemeris_downsampler import build_ephemeris_table, parse_crossing_list

import matplotlib


matplotlib.use('QtAgg')

spice_client = ClientSPICE()

spice_client.KERNEL_LOCATIONS.update(
    {
        "MESSENGER Frames (tf)": {
            "BASE": "https://naif.jpl.nasa.gov/pub/naif/",
            "DIRECTORY": "pds/data/mess-e_v_h-spice-6-v1.0/messsp_1000/data/fk/",
            "PATTERNS": ["msgr_dyn_v600.tf"],
        },
        "MESSENGER": {
            "BASE": "https://naif.jpl.nasa.gov/pub/naif/",
            "DIRECTORY": "pds/data/mess-e_v_h-spice-6-v1.0/messsp_1000/data/spk/",
            "PATTERNS": ["msgr_??????_??????_??????_od431sc_2.bsp"],
        },
    }
)

with spice_client.KernelPool():
    build_ephemeris_table(force_rebuild=True, frame="All")
