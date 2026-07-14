import astropy.units as u
import numpy as np

from matplotlib.lines import Line2D
from hermpy.plotting import Panel
from hermpy.net import ClientSPICE
from hermpy.utils import Constants as c
from astropy.time import Time
from astropy.table import QTable, vstack
from typing import Literal, get_args
import os
import warnings 

from hermpymod.functions.ephemeris_downsampler import parse_crossing_list
from hermpymod.functions.downsampled_positional_data import parse_spice_downsampled
from hermpymod.functions.boundary_models_mod import plot_magnetospheric_boundaries
from hermpymod.functions.encounters import parse_encounters_list


# Define custom Mercury radii unit
mercury_rad = c.MERCURY_RADIUS.to("km")

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

crossing_list = parse_crossing_list()
crossing_times = Time(crossing_list["UTC"]).to_datetime()

encounters_data = parse_encounters_list()
encounter_times = Time(encounters_data["Time Start"]).to_datetime()


"""
Class to plot downsampled positional data, for 2D planar plots
"""

class PlanarplotPanel(Panel):
    def __init__(self, time=[str, str], plane: Literal["X-Y", "X-Z", "Y-Z"] = "X-Y", cylindrical=False, units="Mercury Radii", mercury=True, crossings=False, MP=True, BS=True, encounters=False, add_legend=True, frame = "MSO", alpha=1.0, color='C0'):
        # Initialize the parent Panel class
        super().__init__() 
        
        # Store the plotting configurations as internal state
        self.time = Time(time).to_datetime()
        self.plane = plane
        self.units = units
        self._mercury = mercury
        self._frame = frame
        self._crossings = crossings 
        self._encounters = encounters
        self._legend = add_legend
        self._cylindrical = cylindrical
        self._alpha = alpha
        self._color = color
        self._BS = BS
        self._MP = MP

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
        self._labels = plane_dict[self.plane]
        self._all_labels = [X,Y,Z]

    # Default units are Mercury radii, other options follow astropy units (km, m, etc.)
    def positions(self):

            positions = self.orbit_data[self.poscol].copy()

            for col in positions.keys():
                positions[col] = positions[col].to(self.units)

            return positions
                

    def _plot_on(self, ax):
        # Concentric circles (constant-r gridlines)
        for r in range(0,15):
            theta = np.linspace(0, 2*np.pi, 200)
            ax.plot(r*np.cos(theta), r*np.sin(theta), color='gray', linestyle='--', linewidth=0.5, zorder=0, alpha=0.5)

        # Radial spokes (constant-theta gridlines)  
        for angle_deg in np.arange(0, 360, 30):
            angle = np.radians(angle_deg)
            ax.plot([0, 15*np.cos(angle)], [0, 15*np.sin(angle)], color='gray', linestyle='--', linewidth=0.5, zorder=0)

        if self._cylindrical:
            xlab, ylab, zlab = self._all_labels
            rho = np.sqrt(self.positions()[ylab]**2 + self.positions()[zlab]**2)
            x_table = self.positions()[xlab]
            ax.scatter(x_table, rho, s=0.1, alpha=self._alpha, color=self._color)
            ax.set_xlabel(xlab + r" (" + self.units + r")")
            ax.set_ylabel(f"$\\rho = \\sqrt{{{ylab}^2 + {zlab}^2}}$" + r" (" + self.units + r")")

        else:
            # Try plot from dictionary, if no item exits, return an error
            try:
               xlab, ylab = self._labels
               ax.scatter(self.positions()[xlab], self.positions()[ylab], s=0.1, alpha=self._alpha, color=self._color)
               ax.set_xlabel(xlab + r" (" + self.units + r")")
               ax.set_ylabel(ylab + r" (" + self.units + r")")

            except KeyError:
                raise ValueError("plane must be passed as a string in alphebetical order (X-Y, Y-Z or X-Z)")

        if self._mercury:
            # X^2 + Y^2 + Z^2 = 1 => X^2 + rho^2 = 1
            R = mercury_rad.to(self.units)
            if self._cylindrical:
                t = np.linspace(0, np.pi, 100)
            else:
                t = np.linspace(0, 2* np.pi, 100)
            ax.plot(R*np.cos(t), R*np.sin(t), lw = 2, color='k', label='Mercury')
            ax.set_aspect('equal')
        
        if self._encounters:

            encounter_time_mask = (encounter_times >= self.time_range[0]) & (encounter_times <= self.time_range[-1])
            orbit_encounters = [[i["Time Start"], i["Time End"]] for i in encounters_data[encounter_time_mask]]

            encounters_positions = []

            for encounter in orbit_encounters:
                encounters_positions.append(parse_spice_downsampled(encounter))

            encounters_positions_full = vstack(encounters_positions)

            if self._cylindrical:
                xlab, ylab, zlab = self._all_labels
                encounters_rho = np.sqrt(encounters_positions_full[ylab]**2 + encounters_positions_full[zlab]**2)

                ax.scatter(encounters_positions_full[xlab], encounters_rho, s=0.8, label=f"Encounters: Number of encounters = {len(encounters_positions)}", color="orange", zorder=2.5)
            else:
                xlab, ylab = self._labels
                ax.scatter(encounters_positions_full[xlab], encounters_positions_full[ylab], s=0.8, label=f"Encounters: Number of encounters = {len(encounters_positions)}", color="orange", zorder=2.5)


        if self._crossings:

            # Take crossings that match time range given
            mask = match_times(self.time, crossing_times)

            crossing_list_panel = crossing_list[mask]

            # Dictionary of masks for each crossing type
            mask_dict = {
                "BS_OUT":{"mask": crossing_list_panel["Label"] == "BS_OUT", "color":'yellow'},
                "BS_IN" :{ "mask": crossing_list_panel["Label"] == "BS_IN", "color": 'red'},
                "MP_OUT":{ "mask": crossing_list_panel["Label"] == "MP_OUT","color": 'purple'},
                "MP_IN" :{ "mask": crossing_list_panel["Label"] == "MP_IN", "color": 'blue'},
                    }


            if not self._MP:
                del mask_dict["MP_IN"]
                del mask_dict["MP_OUT"]


            if not self._BS:
                del mask_dict["BS_IN"]
                del mask_dict["BS_OUT"]


            crossing_positions = crossing_list_panel['|R|', 'X MSO', 'Y MSO', 'Z MSO']

            for col in crossing_positions.keys():
                crossing_positions[col] = crossing_positions[col].to(self.units)

            if self._cylindrical:
                xlab, ylab, zlab = self._all_labels
                for mask in mask_dict.keys():
                    rho_crossing = np.sqrt(crossing_positions[ylab][mask_dict[mask]["mask"]]**2 + crossing_positions[zlab][mask_dict[mask]["mask"]]**2)
                    crossing_x = crossing_positions[xlab][mask_dict[mask]["mask"]]
                    ax.scatter(crossing_x, rho_crossing, label=mask, marker='x', s=50, color=mask_dict[mask]["color"])

            else:
                xlab, ylab = self._labels
                for mask in mask_dict.keys():
                    ax.scatter(crossing_positions[xlab][mask_dict[mask]["mask"]],crossing_positions[ylab][mask_dict[mask]["mask"]], label=mask, marker='x', s=50, color=mask_dict[mask]["color"])

            if self.units != "Mercury Radii":
                warnings.warn(
            "Warning, to plot magnetic boundaries you must use units of 'Mercury Radii'",
            UserWarning,
            stacklevel=2,
        )

            plane = {
                    "X-Y": "xy",
                    "Y-Z": "yz",
                    "X-Z": "xz",
                    }
            R = mercury_rad.to(self.units).value

            plot_magnetospheric_boundaries(ax, frame=self._frame, plane=plane[self.plane], add_legend=True, cylindrical=self._cylindrical)

            ax.set_xlim(-10*R, 10*R)
            ax.set_ylim(-10*R, 10*R)

            if self._cylindrical:
                ax.set_ylim(0, 10*R)

            if self._legend:
                ax.legend(loc='best', bbox_to_anchor=(1.2,1.2))


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
