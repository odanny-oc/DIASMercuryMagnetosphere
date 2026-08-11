# Data Handling and Preparation Directory

## Orbit Finder

*Note* it is unnecessary to run any file in this directory for the analysis of crossings and encounters.

### downsample_spice_data.py
Uses the `build_ephemeris_table()` function (from `hermpymod.functions.ephemeris_downsampler`) to download the MESSENGER SPICE kernels. This is done using the custom MESSENGER trajectory and frame kernel patterns with `hermpy`'s `ClientSPICE`. It then downsamples the resulting ephemeris to one point per minute over the entire mission. Since orbit and boundary analysis operates on timescales of hours, this coarse resolution is sufficient and keeps the dataset manageable; the table is cached to disk and only rebuilt when `force_rebuild=True` is passed. When needed by other files, this functions is called in live.

### orbit_periapsis_apoapsis_finder.py
Loads the downsampled ephemeris (`parse_spice_downsampled()`) and the Hollman crossing list (`parse_crossing_list()`), then uses `scipy.signal.find_peaks` on the (negative of the) radial distance `|R|` to locate periapsides and apoapsides, which by definition mark the start/end of each orbit. Periapsis and apoapsis times and positions are rounded and written out as `periapsis_data.csv` and `apoapsis_data.csv`. The file also produces a diagnostic plot of MESSENGER's distance from Mercury with periapsis, apoapsis, and all four crossing types marked, plus 2D planar trajectory plots (X-Y, X-Z, Y-Z) built with the `PlanarplotPanel` class. These datasets can be found on Zenodo (https://zenodo.org/records/21393912).

### orbit_times_analysis.py
Validates the periapsis-defined orbit list by checking that consecutive periapsis-to-periapsis intervals are physically reasonable orbit durations. It flags orbits outside the 7-13 hour range and specifically identifies the eleven ~9-10 hour transition orbits that occur as MESSENGER moves from its initial 12-hour to its later 8-hour orbital period. It plots a log-scaled histogram of orbit lengths (saved as `hollman_orbit_times.svg`), which should show two peaks — one per orbital period — either side of the transition orbits, and plots the radial distance over the transition period itself.

### orbit_data_plotter.py
Groups the Hollman encounter list into orbits using the periapsis times as orbit boundaries, then for whichever orbit number(s) are listed in `plots`, plots the MESSENGER magnetometer data (`mag_data_plotter`) and the full ephemeris/trajectory (`plot_all_ephemeris`) for that orbit, with the orbit window highlighted. Useful for visually inspecting a specific orbit's crossings and encounters against the underlying MAG data.

## List of Magnetospheric Boundary Encounters

An encounter is defined as the time spent in the vicinity of a magnetospheric boundary (bow shock/magnetopause). This directory prepares the Sun et al. and Philpott et al. encounter lists into a common four-column format (Time Start, Time End, Label, Encounter Duration) so they can be directly compared against the Hollman list elsewhere in the repository.

### philpott_data.py
Reads the Philpott et al. `supporting_table_S1.tab` boundary-crossing table(https://borealisdata.ca/dataset.xhtml?persistentId=doi:10.5683/SP2/1U6FEO), drops the two data-gap boundary numbers (9 and 10), and relabels the eight numeric `Boundarynumber` categories into the (BS_IN outer/inner, MP_IN inner/outer, MP_OUT inner/outer, BS_OUT inner/outer) scheme. It writes a flat crossing list (`philpott_crossings_list_2020.csv`), then pairs up the inner/outer boundary crossings into encounters and writes the result to `philpott_encounter_list_2020.csv` in the standard four-column encounter format.

### sun_data.py
Downloads the four Sun et al. bow shock/magnetopause in/out `.txt` files from Zenodo via `curl` (https://zenodo.org/records/18236915), parses each fixed-width line with `parse_bsi_line()` (start time, end time, label, qualifier) using the `sun_data_parser()` function, computes encounter duration, and concatenates and time-orders all four files into a single `sun_2023_crossing.csv`, again in the standard four-column encounter format.

Either of the directories, "encounters" or "crossings", can now be run depending on the desire of the user.
