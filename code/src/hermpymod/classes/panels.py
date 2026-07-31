import astropy.units as u
import numpy as np

from matplotlib.lines import Line2D
from hermpy.plotting import Panel
from hermpy.net import ClientSPICE
from hermpy.utils import Constants as c
from astropy.time import Time
from astropy.table import QTable, vstack, hstack, Table
from typing import Literal, get_args
import os
import warnings 

from hermpymod.functions.ephemeris_downsampler import parse_crossing_list, parse_spice, time_array, abs_r
from hermpymod.functions.downsampled_positional_data import parse_spice_downsampled
from hermpymod.functions.boundary_models_mod import plot_magnetospheric_boundaries
from hermpymod.functions.encounters import parse_encounters_list


# Define custom Mercury radii unit
mercury_rad = c.MERCURY_RADIUS.to("km")

Zd = c.DIPOLE_OFFSET.to("km")


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


def polar_grid(ax):
    # Concentric circles (constant-r gridlines)
    for r in range(0,15):
        theta = np.linspace(0, 2*np.pi, 200)
        ax.plot(r*np.cos(theta), r*np.sin(theta), color='gray', linestyle='--', linewidth=0.5, zorder=0, alpha=0.2)

    # Radial spokes (constant-theta gridlines)  
    for angle_deg in np.arange(0, 360, 30):
        angle = np.radians(angle_deg)
        ax.plot([0, 15*np.cos(angle)], [0, 15*np.sin(angle)], color='gray', linestyle='--', linewidth=0.5, zorder=0, alpha=0.2)

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
    def __init__(self, time=[str, str], plane: Literal["X-Y", "X-Z", "Y-Z", "All"] = "X-Y", cylindrical=False, units="Mercury Radii", mercury=True, crossings=None, MP=True, BS=True, encounters=None, add_legend=True, frame = "MSO", alpha=1.0, color='C0', scatter=True, label=None, grid=False, downsampled=True, resolution=None):
        # Initialize the parent Panel class
        super().__init__() 
        
        # Store the plotting configurations as internal state
        self.time = Time(time).to_datetime()
        self.plane = plane
        self.units = units
        self._mercury = mercury
        self._frame = frame
        self._grid = grid
        self._crossings = crossings
        self._encounters = encounters
        self._legend = add_legend
        self._cylindrical = cylindrical
        self._alpha = alpha
        self._color = color
        self._BS = BS
        self._MP = MP
        self._scatter = scatter
        self._label = label

        if downsampled:
            orbit_data = parse_spice_downsampled(time_range=time)

        else:

            if resolution==None:
                raise ValueError("'resolution' not defined, select resolution of sampling in seconds.")
            else:
                time_range = time_array(Time(time[0]).to_datetime(), Time(time[-1]).to_datetime(), resolution)

            orbit_data = parse_spice(time_range, units=self.units, frame=self._frame)

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
                "All": [X, Y ,Z],
                }

        self._labels = plane_dict[self.plane]
        self._all_labels = [X,Y,Z]

    # Default units are Mercury radii, other options follow astropy units (km, m, etc.)
    def positions(self):

            positions = self.orbit_data[self.poscol].copy()

            for col in positions.keys():
                positions[col] = positions[col].to(self.units)

            return positions
                
    def plot_trajectories(self , ax, labels):
    
        def traj_plotter(x,y):
            if self._scatter:
                return ax.scatter(x, y, s=0.1, alpha=self._alpha, color=self._color, label=self._label)
            else:
                return ax.plot(x, y, alpha=self._alpha, color=self._color, label = self._label)

        # Try plot from dictionary, if no item exits, return an error
        if self._cylindrical:
            xlab, ylab, zlab = self._all_labels
            rho = np.sqrt(self.positions()[ylab]**2 + self.positions()[zlab]**2)
            x_table = self.positions()[xlab]
            traj_plotter(x_table, rho)
            ax.set_xlabel(xlab + r" (" + self.units + r")")
            ax.set_ylabel(f"$\\rho = \\sqrt{{{ylab}^2 + {zlab}^2}}$" + r" (" + self.units + r")")
        else:
            try:
               xlab, ylab = labels
               traj_plotter(self.positions()[xlab], self.positions()[ylab])
               ax.set_xlabel(xlab + r" (" + self.units + r")")
               ax.set_ylabel(ylab + r" (" + self.units + r")")

            except KeyError:
                raise ValueError("plane must be passed as a string in alphebetical order (X-Y, Y-Z or X-Z)")


    def plot_mercury(self, ax, label):
        # X^2 + Y^2 + Z^2 = 1 => X^2 + rho^2 = 1
        R = mercury_rad.to(self.units)
        t = np.linspace(0, 2* np.pi, 100)
        t2 = np.linspace(0, np.pi, 100)

        if self._frame == "MSO":
            if self._cylindrical:
                y = abs(R*np.sin(t))
            else:
                y = R*np.sin(t)

            ax.plot(R*np.cos(t), y, lw = 2, color='k', label='Mercury', zorder=1)
            ax.set_aspect('equal')

        elif self._frame == "MSM":
            if self._cylindrical:
                y = abs(R*np.sin(t))
                ax.plot(R*np.cos(t), y + Zd.to(self.units), lw = 2, color='k', label='Mercury', zorder=1)

            elif "Z" in label[-1]:
                y = R*np.sin(t) + Zd.to(self.units)

            else:
                y = R*np.sin(t)

            if "Z" in label[-1]:
                y = R*np.sin(t) + Zd.to(self.units)

            ax.plot(R*np.cos(t), y, lw = 2, color='k', label='Mercury', zorder=1)
            ax.set_aspect('equal')



    def plot_encounters(self, ax, labels):

        encounter_start_times = Time(self._encounters["Time Start"]).to_datetime()
        encounter_end_times = Time(self._encounters["Time End"]).to_datetime()
        encounter_time_start_mask = (encounter_start_times >= self.time[0]) & (encounter_start_times <= self.time[-1])
        encounter_time_end_mask = (encounter_end_times >= self.time[0]) & (encounter_end_times <= self.time[-1])

        encounter_time_mask = encounter_time_start_mask | encounter_time_end_mask

        orbit_encounters = [[i["Time Start"], i["Time End"]] for i in self._encounters[encounter_time_mask]]

        encounters_positions = []

        for encounter in orbit_encounters:
            encounters_positions.append(parse_spice_downsampled(encounter))

        encounters_positions_full = vstack(encounters_positions)
        self.encounters = encounters_positions_full

        if self._cylindrical:
            xlab, ylab, zlab = self._all_labels
            encounters_rho = np.sqrt(encounters_positions_full[ylab]**2 + encounters_positions_full[zlab]**2)

            ax.scatter(encounters_positions_full[xlab], encounters_rho, s=0.8, label=f"Encounters: Number of encounters = {len(encounters_positions)}", color="orange", zorder=4)
        else:
            xlab, ylab = labels
            ax.scatter(encounters_positions_full[xlab], encounters_positions_full[ylab], s=0.8, label=f"Encounters: Number of encounters = {len(encounters_positions)}", color="orange", zorder=4)


    def plot_crossings(self, ax, labels):

        crossing_times = Time(self._crossings["UTC"]).to_datetime()
        # Take crossings that match time range given
        mask = match_times(self.time, crossing_times)

        crossing_list_panel = self._crossings[mask]

        self.crossings = crossing_list_panel


        # Dictionary of masks for each crossing type
        mask_dict = {
            f"BS_OUT\nNumber of crossings {len(crossing_list_panel)}":{"mask": (crossing_list_panel["Label"] == "BS_OUT") | (crossing_list_panel["Label"] == "BSO") , "color":'yellow'},
            "BS_IN" :{ "mask": (crossing_list_panel["Label"] == "BS_IN") | (crossing_list_panel["Label"] == "BSI"), "color": 'red'},
            "MP_OUT":{ "mask": (crossing_list_panel["Label"] == "MP_OUT") | (crossing_list_panel["Label"] == "MPO"), "color": 'purple'},
            "MP_IN" :{ "mask": (crossing_list_panel["Label"] == "MP_IN") | (crossing_list_panel["Label"] == "MPI"), "color": 'blue'},
                }


        if not self._MP:
            MP_keys = [i for i in mask_dict.keys() if "MP" in i]
            for key in MP_keys:
                del mask_dict[key]


        if not self._BS:
            BS_keys = [i for i in mask_dict.keys() if "BS" in i]
            for key in BS_keys:
                del mask_dict[key]


        if self._frame == "MSM":
            try:
                distance = [abs_r([x,y,z]) for x,y,z in zip(crossing_list_panel["X MSM"], crossing_list_panel["Y MSM"] , crossing_list_panel["Z MSM"] + Zd.to(self.units).value)]
                crossing_positions = QTable({
                    "|R|" : distance,
                    "X MSO" : crossing_list_panel["X MSM"],
                    "Y MSO" : crossing_list_panel["Y MSM"],
                    "Z MSO" : crossing_list_panel["Z MSM"]
                    })

            except KeyError:
                distance = [abs_r([x,y,z]) for x,y,z in zip(crossing_list_panel["X MSO"], crossing_list_panel["Y MSO"] , crossing_list_panel["Z MSO"] + Zd.to(self.units).value)]
                
                crossing_positions = QTable({
                    "|R|" : distance,
                    "X MSO" : crossing_list_panel["X MSO"],
                    "Y MSO" : crossing_list_panel["Y MSO"],
                    "Z MSO" : crossing_list_panel["Z MSO"] - Zd.to(self.units).value})

        elif self._frame == "MSO":
            crossing_positions = crossing_list_panel['|R|', 'X MSO', 'Y MSO', 'Z MSO']

        # for col in crossing_positions.keys():
        #     crossing_positions[col] = crossing_positions[col].to(self.units)

        if self._cylindrical:
            xlab, ylab, zlab = self._all_labels
            for mask in mask_dict.keys():
                rho_crossing = np.sqrt(crossing_positions[ylab][mask_dict[mask]["mask"]]**2 + crossing_positions[zlab][mask_dict[mask]["mask"]]**2)
                crossing_x = crossing_positions[xlab][mask_dict[mask]["mask"]]
                ax.scatter(crossing_x, rho_crossing, label=mask, marker='x', s=50, color=mask_dict[mask]["color"], zorder = 5)

        else:
            xlab, ylab = labels
            for mask in mask_dict.keys():
                ax.scatter(crossing_positions[xlab][mask_dict[mask]["mask"]],crossing_positions[ylab][mask_dict[mask]["mask"]], label=mask, marker='x', s=50, color=mask_dict[mask]["color"], zorder=5)

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

        plot_magnetospheric_boundaries(ax, frame=self._frame, plane=plane[self.plane], add_legend=True, cylindrical=self._cylindrical, zorder=3)

        ax.set_xlim(-10*R, 10*R)
        ax.set_ylim(-10*R, 10*R)

        if self._cylindrical:
            ax.set_ylim(0, 10*R)

        if self._legend:
            ax.legend(loc='best', bbox_to_anchor=(1.2,1.2))

    def _plot_on(self, ax):
        # Plot grid lines
        if self.plane == "All":
            if len(ax) != 4:
                raise ValueError("Not enough axis to plot 'All'")
            else:
                plot_config = [
                        (ax[0], "X-Y", [self._labels[0], self._labels[1]], False),
                        (ax[1], "X-Z", [self._labels[0], self._labels[2]], False),
                        (ax[2], "Y-Z", [self._labels[1], self._labels[2]], False),
                        (ax[3], "X-Y", [self._labels[0], self._labels[1]], True),
                        ]
                for axis, plane, label, cylindrical in plot_config:
                    if not cylindrical and self._grid:
                        polar_grid(axis)
                        axis.set_xlim(-5,5)
                        axis.set_ylim(-5,5)
                    self._cylindrical = cylindrical
                    self.plane = plane
                    self.plot_trajectories(axis, label)
                    if self._mercury:
                        self.plot_mercury(axis, label)
                    if self._crossings != None:
                        self.plot_crossings(axis, label)
                    if self._encounters != None:
                        self.plot_encounters(axis, label)
                
        else:
            if self._grid:
                polar_grid(ax)
            self.plot_trajectories(ax, labels=self._labels)

            if self._mercury:
                self.plot_mercury(ax, labels=self._labels)
            
            if self._encounters != None:
                self.plot_encounters(ax, labels=self._labels)

            if self._crossings != None:
                self.plot_crossings(ax, labels=self._labels)


class HistogramPanel(Panel):
        def __init__(self, data, bins='auto', color='C0', minmax=False, average=True, label=None, zorder=2):
            # 1. Properly initialize the parent Panel class
            super().__init__() 
            self._data = data
            self._bins = bins
            self._color = color
            self.average = np.average(data)
            self._average = average
            self._hist, _ = np.histogram(data, bins=bins)
            self._minmax = minmax
            self._label = label
            self._zorder = zorder

        def _plot_on(self, ax):
            ax.hist(self._data, bins=self._bins, edgecolor='k', color=self._color, label=self._label, zorder=self._zorder)
            if self._average:
                ax.axvline(self.average, ls='--', color='r', label=f'Average = {self.average:.2f}')
            plot_handles, _ = ax.get_legend_handles_labels()
            handles = [Line2D([], [], color = 'none', label = f"Number of points = {len(self._data)} \nMedian {np.median(self._data):.2f}\nStandard Deviation {np.std(self._data):.2f}")]
            if self._minmax:
                minmax = [Line2D([], [], color = 'none', label = f"Max data = {np.max(self._data):.2f} \nMin data {np.min(self._data):.2f}")]
                handles += minmax
            ax.legend(handles= plot_handles + handles)
