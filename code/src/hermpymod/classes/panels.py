import astropy.units as u
import numpy as np
from hermpy.plotting import Panel
from hermpy.net import ClientSPICE
from astropy.time import Time
import spiceypy as spice
from hermpy.utils import Constants as c
import datetime as dt
from matplotlib.lines import Line2D
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from functions.ephemeris_downsampler import parse_spice_downsampled, parse_spice, parse_crossing_list, abs_r


# Define custom Mercury radii unit
mercury_rad = c.MERCURY_RADIUS.to("km")
R_M = u.def_unit("R_M", mercury_rad)
u.add_enabled_units(R_M)

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


"""
Function to create masks for given time intervals
"""
def match_times(time_array, compare_array):
    mask = (time_array[0] <= compare_array) & ((time_array[-1] >= compare_array))
    return mask


home_dir = os.getenv('HOME')
data_dir = os.path.join(home_dir, '.ephemeris_data/')
os.makedirs(data_dir, exist_ok = True)


"""
Class to plot downsampled positional data, for 2D planar plots
"""

class PlanarplotPanel(Panel):
    def __init__(self, time=[str, str], plane=None, units="R_M", mercury=True, crossings=False):
        # Initialize the parent Panel class
        super().__init__() 
        
        # Store the plotting configurations as internal state
        self.time = Time(time).to_datetime()
        self.plane = plane
        self.units = units
        self._mercury = mercury
        self._crossings = crossings 

        orbit_data = parse_spice_downsampled(time_range=time)

        self.orbit_data = orbit_data
        self.time_range = self.orbit_data["UTC"].to_datetime()
        self.poscol = [i for i in self.orbit_data.keys() if "X" in i or "Y" in i or "Z" in i]

        # Dictionary to plot based on given plane
        X = [i for i in self.poscol if "X" in i][0]
        Y = [i for i in self.poscol if "Y" in i][0]
        Z = [i for i in self.poscol if "Z" in i][0]

        plane_dict = {
                "X-Y": [X, Y],
                "Y-Z": [Y, Z],
                "X-Z": [X, Z],
                }
        self.labels = plane_dict[self.plane]

    # Default units are Mercury radii, other options follow astropy units (km, m, etc.)
    def positions(self):

            positions = self.orbit_data[self.poscol].copy()

            for col in positions.keys():
                positions[col] = positions[col].to(self.units)

            return positions
                

    def _plot_on(self, ax):
        # Try plot from dictionary, if no item exits, return an error
        try:
           xlab, ylab = self.labels
           ax.scatter(self.positions()[xlab], self.positions()[ylab], s=0.1)
           ax.set_xlabel(xlab + r" (" + self.units + r")")
           ax.set_ylabel(ylab + r" (" + self.units + r")")

        except KeyError:
            raise ValueError("plane must be passed as a string in alphebetical order (X-Y, Y-Z or X-Z)")

        if self._mercury:
            R = mercury_rad.to(self.units)
            t = np.linspace(0, 2 * np.pi, 100)
            ax.plot(R*np.cos(t), R*np.sin(t), lw = 2, color='k', label='Mercury')
            ax.set_aspect('equal')

        if self._crossings:

            crossing_list = parse_crossing_list()
            
            crossing_times = Time(crossing_list["UTC"]).to_datetime()

            # Take crossings that match time range given
            mask = match_times(self.time, crossing_times)

            crossing_list = crossing_list[mask]
            crossing_times = crossing_times[mask]

            xlab, ylab = self.labels

            # Dictionary of masks for each crossing type
            mask_dict = {
                "BS_OUT":{"mask": crossing_list["Label"] == "BS_OUT", "color":'yellow'},
                "BS_IN" :{ "mask": crossing_list["Label"] == "BS_IN", "color": 'red'},
                "MP_OUT":{ "mask": crossing_list["Label"] == "MP_OUT","color": 'purple'},
                "MP_IN" :{ "mask": crossing_list["Label"] == "MP_IN", "color": 'blue'},
                    }

            crossing_positions = crossing_list['|R|', 'X MSO', 'Y MSO', 'Z MSO']

            for col in crossing_positions.keys():
                crossing_positions[col] = crossing_positions[col].to(self.units)

            for mask in mask_dict.keys():
                ax.scatter(crossing_positions[xlab][mask_dict[mask]["mask"]],crossing_positions[ylab][mask_dict[mask]["mask"]], label=mask, marker='x', s=50, color=mask_dict[mask]["color"])
            ax.legend()


class HistogramPanel(Panel):
        def __init__(self, data, bins='auto', color='C0', minmax=False):
            # 1. Properly initialize the parent Panel class
            super().__init__() 
            self._data = data
            self._bins = bins
            self._color = color
            self.average = np.average(data)
            self._hist, _ = np.histogram(data, bins=bins)
            self._minmax = minmax

        def _plot_on(self, ax):
            ax.hist(self._data, bins=self._bins, edgecolor='k', color=self._color)
            ax.axvline(self.average, ls='--', color='r', label=f'Average = {self.average:.2f}')
            plot_handles, _ = ax.get_legend_handles_labels()
            handles = [Line2D([], [], color = 'none', label = f"Number of points = {len(self._data)} \nMedian {np.median(self._data):.2f}\nStandard Deviation {np.std(self._data):.2f}")]
            if self._minmax:
                minmax = [Line2D([], [], color = 'none', label = f"Max data = {np.max(self._data):.2f} \nMin data {np.min(self._data):.2f}")]
                handles += minmax
            ax.legend(handles= plot_handles + handles)
