# Crossings Directory

A crossing is recorded as the time when the spacecraft crosses the bow shock or magnetopause, as classified by the Hollman 2025 region classifier. Files in this directory can be run in any order.

## Generic Analysis & Visualisations (./01_generic_analysis_and_visualisations/)

### all_bs_crossings.py / all_mp_crossings.py
Plot every bow shock (or magnetopause) crossing in the Hollman list in the X-ρ (cylindrical MSO) plane, coloured by azimuthal angle φ. For each crossing, `encounter_finder()` groups it with its encounter, and an error bar shows how far the spacecraft moved (in X and ρ) over the course of that encounter. This gives a sense of how much a single "crossing" point actually smears across the boundary. `plot_magnetospheric_boundaries()` overlays the average model boundary from Winslow et al. (https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/jgra.50237) for reference.

### dt_between_crossings.py
Computes the time between consecutive crossings via `delta_crossing_times()`, split into all crossings, bow-shock-only, and magnetopause-only. Plots log-scaled histograms of these time gaps; short gaps flag periods of rapid, highly variable boundary motion. The short gaps are expected to be associated withe the variability during encounters, while the longer gaps could be associated with variability between encounters.

### per_orbit_by_type.py
Uses `orbit_data()` to group crossings by orbit, then for each of the four crossing types (BS, MP, IN, OUT) counts how many occur per orbit and histograms the result. Includes a sanity check that the total crossings counted across all orbits matches the length of the full crossing list.

### total_per_orbit.py
Same orbit grouping, but histograms the total number of crossings (of any type) per orbit on a log scale, giving an overall picture of how "busy" a typical orbit is in terms of boundary crossings.

## Extrema & Outlier Analysis (./02_extrema_and_outlier_analysis/)

### most_crossings_per_orbit.py
Finds the orbit(s) with the highest number of bow shock crossings, then plots the full ephemeris/trajectory and MAG data for the orbit with the most crossings and the orbit with the second-most, with the encounter and orbit windows highlighted, so the underlying magnetic field behaviour driving the high crossing count can be inspected directly.

### odd_crossing.py
A well-formed orbit should contain an even number of BS and MP crossings. This file finds every orbit where that isn't true, and plots the MAG data for each so the odd count can be explained. We find that these are due to either a data gap or a misplaced/missing crossing in the Hollman list.

### unphysical.py
Filters for crossings whose label contains the Hollman list's own `UNPHYSICAL` flag, identifies which orbits they belong to, and plots the MAG data around each of the (five) unphysical crossings for manual inspection. The `UNPHYSICAL` flag is a cconsequence of when the region classifier identifies two touching regions which are unphysical i.e solar wind -> magnetosphere, and magnetosphere -> solar wind.

## Crossings Close to Apoapsis & Periapsis (./03_crossings_close_to_orbital_extrema/)

### dt_between_type_apoapsis_periapsis.py
For every bow shock crossing, finds the time to the nearest apoapsis and, separately, to the next apoapsis ahead of it; does the equivalent for magnetopause crossings against periapsis. Results are split by the mission's 12-hour and 8-hour orbit epochs (orbit index 800 marks the transition). Histograms these time differences, and for crossings that sit close to apoapsis (bow shock) or unusually close to periapsis (magnetopause), plots the trajectory around that crossing together with a fitted boundary model (`boundary_fitter`) to visualise how compressed the boundary appears to be at that point.

### mp_close_to_periapsis.py
A focused, single-example version of the above: plots the full ephemeris/trajectory and MAG data for a specific hard-coded time window (30 November 2013) in which a magnetopause crossing occurs very close to periapsis.

## Magnetosheath Traversals (./04_magnetosheath_traversals/)

A magnetosheath traversal is defined from the first and last crossings bordering the magnetosheath. The `delta_t_magenetosheath()` function  that walks the Hollman crossing list and pairs up MP_OUT→BS_OUT crossings as "outbound" magnetosheath traversals and BS_IN→MP_IN crossings as "inbound" traversals.

### time_in_magnetosheath.py
Computes the duration of every inbound/outbound traversal, plots the trajectory of the 100 shortest traversals of each type, and histograms all traversal durations to estimate the total time MESSENGER spent in the magnetosheath.

<img src="../../plots_and_images/inbound_magnetosheath_traversals.svg" width="250">
<img src="../../plots_and_images/outbound_magnetosheath_traversals.svg" width="250">


### shortest_traversals_hist.py
Takes the 100 shortest traversals of each type and spatially bins the spacecraft's position along each one (X-Y, X-Z, Y-Z and cylindrical MSM planes) to produce 2D histograms of the *average traversal duration* per spatial bin, revealing where in space the shortest (most compressed) traversals tend to occur.

### traversals_per_hour_hist.py
Similar spatial binning, but computes the *rate* of traversals per hour of MESSENGER residence time in each spatial bin (using the downsampled ephemeris for the residence-time normalisation), for both inbound and outbound traversals.

### radial_distance_to_boundary.py
For the shortest 100 traversals, computes each crossing's radial distance from Mercury and compares it to the radially-projected average boundary distance at that point (`delta_r`). All distances are measured through the angle $0 \leq \theta=\arctan{\rho/X} \leq 180]$ of the location of the crossing. We plot the raw distances and a 10-day boxcar-smoothed running average of the difference (the "compression distance") over the mission. Also produces a 2D histogram of _"compression factor"_, defined as 
$$\frac{(\text{radial distance from boundary to MESSENGER})}{\text{radial distance from Mercury to boundary}},$$
against boundary angle $\theta$, coloured by mean traversal duration.

### compression_factor.py
A near-identical companion to `radial_distance_to_boundary.py`, focused specifically on the dimensionless compression factor (`delta_r / boundary_distance`, rather than the raw distance), plotted as a 10-day boxcar-smoothed time series for each boundary/direction combination.

### grazing_angle_compression_factor.py
Extends the compression factor calculation by also computing the grazing angle (`get_grazing_angle()`, the angle between the boundary normal and MESSENGER's velocity) for each crossing, and produces a 2D histogram of compression factor vs. boundary angle $\theta$, coloured by mean grazing angle.

### traversal_grazing_angle.py
For the 100 shortest inbound and outbound traversals, computes the grazing angle at the bow shock and magnetopause crossings (`get_grazing_angle(..., return_vectors=True)`), then plots the trajectories coloured by grazing angle with the boundary normal and MESSENGER velocity vectors overplotted as quivers at each crossing, for direct visual interpretation of how perpendicular the spacecraft's approach to the boundary was.

### crossing_varibility.py
Looks at the 50 shortest traversals in detail: plots their trajectories in MSM coordinates alongside their neighbouring (adjacent-orbit) crossings for comparison, then fits a boundary model (`boundary_fitter`) to the first five traversal/adjacent-crossing pairs and overlays the fitted magnetopause or bow shock surface on the trajectory plot.

### mag_data_magneosheath_traversals.py
Plots the MESSENGER magnetometer data across the shortest magnetosheath traversals (zoomed to the traversal itself, with a 3-hour context window either side) with the crossings marked, to visually inspect the magnetic field signature during the shortest, most compressed sheath passages.

### spatial_mag_data.py
Plots the trajectory of every 10th traversal among the shortest 100 (for both directions) in the X-Y, X-Z, Y-Z and cylindrical MSM planes, with the trajectory colour-coded by total magnetic field strength |B|, giving a combined spatial/field-strength view of the shortest traversals.
