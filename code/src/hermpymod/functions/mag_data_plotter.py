import numpy as np

from astropy.table import QTable
from astropy.time import Time
from sunpy.time import TimeRange

from hermpy.plotting import  TimeseriesPanel
from hermpy.data import parse_messenger_mag
from hermpy.net import ClientMESSENGER
from hermpy.utils import Constants as c

from hermpymod.functions.ephemeris_downsampler import parse_crossing_list
from hermpymod.functions.crossings_plotter import plot_crossings, plot_encounters


c = ClientMESSENGER()


def total_mag_field(data):
    total_data = data["UTC", "Bx", "By", "Bz"]
    total_mag_data = np.sqrt(total_data["Bx"]**2 + total_data["By"]**2 + total_data["Bz"]**2)
    mag_field_data = QTable(data=[data["UTC"], total_mag_data], names=["UTC", "|B|"])
    return total_data, mag_field_data


def mag_data_plotter(time, crossings=True, encounters=True):

    if not isinstance(time[0], str):
        t_start = time[0]
        t_end = time[-1]
        str_times = [t_start.isoformat(), t_end.isoformat()]
    else:
        str_times = [time[0], time[-1]]
        t_start = Time(time[0]).to_datetime()
        t_end =  Time(time[-1]).to_datetime()


    time_range = TimeRange(t_start, t_end)

    c.query(time_range, "MAG")
    mag_data_encounter = c.fetch()

    mag_table : QTable = parse_messenger_mag(mag_data_encounter, time_range)

    directional_mag_data, total_mag_data = total_mag_field(mag_table)
    total_mag_plot = TimeseriesPanel(total_mag_data)
    directional_mag_plot = TimeseriesPanel(directional_mag_data)
    
    mag_plot = directional_mag_plot + total_mag_plot

    fig_mag, ax_mag = mag_plot.plot(show=False)

    fig_mag.suptitle(f"MESSENGER MAG data, taken from {str_times[0][:10]}-{str_times[-1][:10]}")
    lines = ax_mag[1].get_lines()
    lines[0].set_color('k')
    ax_mag[0].axhline(0, ls='--', color='k', label='Zero line')
    ax_mag[0].set_xlabel("Time (UTC)")

    if crossings:
        plot_crossings(t_start, t_end, ax_mag)

    if encounters:
        plot_encounters(t_start, t_end, ax_mag)

    for ax in ax_mag:
        ax.legend()

    return fig_mag, ax_mag
