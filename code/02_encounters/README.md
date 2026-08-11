# Encounters Directory

An encounter is a collection of crossings of the same type and direction (e.g. first bow shock in to last bow shock in), as defined using the trajectory direction field in the Hollman list.

## Generic Analysis & Visualisation (./01_generic_analysis_and_visualisation/)

### encounters_per_orbit.py
**Note** to run this file one needs to have ran the scripts in the `../00_data_handling/magnetospheric_boundary_encounters` directory, to properly prepare the datasets for comparison.

Groups the Hollman, Sun et al., and Philpott et al. encounter lists into orbits (using the periapsis-defined orbit boundaries) by matching each encounter's start/end times to an orbit window, then histograms the number of encounters per orbit for each dataset. Also reports how many encounters were "overcounted" i.e assigned to two orbits because they straddle an orbit boundary.

### dt_between_encounters.py
Computes the time between consecutive encounters of the same type (BS-BS or MP-MP, via `dt_encounters()`) for the Hollman, Sun, and Philpott encounter lists, and histograms these gaps alongside the raw encounter durations. Defines a `CDFPanel` class to plot _cumulative distribution functions_ of the BS and MP inter-encounter times for each pair of datasets, annotated with a two-sample Kolmogorov-Smirnov statistic and Wasserstein distance (work), to statistically compare how the three encounter lists differ from one another. A final pass (`mode='diff'`) also looks at the time between *different*-type encounters (BS followed by MP, and vice versa).

### crossings_per_encounter.py
For every encounter, uses `arclength_midpoint()` to find the spacecraft's position at the arc-length midpoint of the encounter (interpolating along the downsampled ephemeris), then scatters these midpoint positions (in MSM coordinates, across X-Y, X-Z, Y-Z and cylindrical planes) coloured by the number of crossings in that encounter, with error bars showing how far the spacecraft moved over the encounter. Also histograms the number of BS, MP, and total crossings per encounter.

## Extrema & Outlier Analysis (./02_extrema_and_outlier_analysis/)

### even_encounters.py
A well-formed encounter must start and end with a crossing of the same label, and so should always contain an odd number of crossings. This finds every encounter with an even crossing count and plots the MAG data (with the encounter window highlighted) for a representative sample of them, to inspect what's driving the anomaly.

### overlap_encounters.py / philpott_overlap_encounters.py
Groups encounters into orbits and checks for cases where the last encounter of one orbit is the same encounter as the first encounter of the next orbit i.e. it straddles the orbit boundary and gets counted for both, which explains the overcounting reported in `encounters_per_orbit.py`. For a sample of these overlapping encounters, plots the MAG data with both the encounter and the two adjacent orbit windows highlighted. `overlap_encounters.py` runs this analysis on the Hollman list; `philpott_overlap_encounters.py` applies the same analysis to the Philpott et al. encounter list.

## Spatial Histograms (./03_spatial_histograms/)

### crossings_per_encounter_2d_hist.py
Spatially bins the downsampled ephemeris (X-Y, X-Z, Y-Z, and cylindrical MSM planes) and, for bow shock and magnetopause encounters separately, computes the average number of crossings-per-encounter within each spatial bin, producing 2D colour maps that show where in space encounters tend to be made up of unusually many (or few) individual crossings.

### encounter_per_hour.py
Uses the same spatial-binning approach as the traversal-rate analysis in `04_magnetosheath_traversals`, but for encounters rather than traversals: bins MESSENGER's total residence time and the number of BS/MP encounters spatially, and divides to get the *rate* of encounters per hour of residence in each bin, plotted alongside a residence-time-only reference histogram.
