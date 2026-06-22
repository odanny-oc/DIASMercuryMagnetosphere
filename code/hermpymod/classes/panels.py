import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import QTable, vstack
from hermpy.plotting import Panel, TimeseriesPanel
from hermpy.net import ClientSPICE
from astropy.time import Time
import spiceypy as spice
from hermpy.utils import Constants as c
import datetime as dt
from matplotlib.lines import Line2D
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from data.ephemeris_downsampler import parse_spice_downsampled

# Implement more robust data handling

# Define custom Mercury radii unit
mercury_rad = c.MERCURY_RADIUS.to("km")
R_M = u.def_unit("R_M", mercury_rad)

u.add_enabled_units(R_M)

"""
Class to plot downsampled positional data, for 2D planar plots
"""

class PlanarplotPanel(Panel):
    def __init__(self, time=[str, str], plane=None, units="R_M", mercury=True):
        # 1. Properly initialize the parent Panel class
        super().__init__() 
        
        # 2. Store the plotting configurations as internal state
        self.time = Time(time).to_datetime()
        self.plane = plane
        self.units = units
        self._mercury = mercury

        orbit_data = parse_spice_downsampled(time_range=time)

        self.orbit_data = orbit_data
        self.time_range = self.orbit_data["UTC"].to_datetime()

    # Default units are Mercury radii, other options follow astropy units (km, m, etc.), vector=False gives vector magnitude
    def positions(self, units = "R_M"):

            positions = self.orbit_data["X MSO", "Y MSO", "Z MSO"]

            positions["X MSO"] = positions["X MSO"].to(units)
            positions["Y MSO"] = positions["Y MSO"].to(units)
            positions["Z MSO"] = positions["Z MSO"].to(units)

            return positions
                

    def _plot_on(self, ax):
        if self.plane == "X-Y":
            ax.scatter(self.positions(units=self.units)["X MSO"], self.positions(units=self.units)["Y MSO"], s=0.1)
            ax.set_xlabel("X MSO " + r"(" + self.units + r")")
            ax.set_ylabel("Y MSO " + r"(" + self.units + r")")

        elif self.plane == "Y-Z":
            ax.scatter(self.positions(units=self.units)["Y MSO"], self.positions(units=self.units)["Z MSO"], s=0.1)
            ax.set_xlabel("Y MSO " + r"(" + self.units + r")")
            ax.set_ylabel("Z MSO " + r"(" + self.units + r")")

        elif self.plane == "X-Z":
            ax.scatter(self.positions(units=self.units)["X MSO"], self.positions(units=self.units)["Z MSO"], s=0.1)
            ax.set_xlabel("X MSO " + r"(" + self.units + r")")
            ax.set_ylabel("Z MSO " + r"(" + self.units + r")")

        else:
            raise ValueError("plane must be passed as a string in alphebetical order (X-Y, Y-Z or X-Z)")
        if self._mercury:
            t = np.linspace(0, 2 * np.pi, 100)
            ax.plot(np.cos(t), np.sin(t), lw = 2, color='k', label='Mercury')


class HistogramPanel(Panel):
        def __init__(self, data, bins='auto', color='C0'):
            # 1. Properly initialize the parent Panel class
            super().__init__() 
            self._data = data
            self._bins = bins
            self._color = color
            self.average = np.average(data)

        def _plot_on(self, ax):
            ax.hist(self._data, bins=self._bins, edgecolor='k', color=self._color)
            ax.axvline(self.average, ls='--', color='r', label=f'Average = {self.average:.2f}')
            plot_handles, _ = ax.get_legend_handles_labels()
            handles = [Line2D([], [], color = 'none', label = f"Number of points = {len(self._data)}")]
            ax.legend(handles= plot_handles + handles)
