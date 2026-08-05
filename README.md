# DIAS Magnetosphere Group Computing Research Internship

This repository details the progress made on the analysis of the magnetic boundary crossings made by the MESSENGER spacecraft from 2011 to 2015, for the DIAS Computing Research Internship.

To start, run the orbit_times_analysis file. This reads the "peaks_data.csv", which is a list of times when MESSENGER was at periapsis. This file verifies that all the orbits are the correct times as well as investigates the outliers of the orbits, namely the eleven ~9.5hr transition orbits that occur due to MESSENGER's two-stage transition from 12hr to 8hr.

Next, one runs the sun_data.py and philpott_list.py files to download the alternative encounter data sets for comparison later. To run philpott_list.py, one must download the Philpott et al. list of encounters https://borealisdata.ca/dataset.xhtml?persistentId=doi:10.5683/SP2/1U6FEO (the supporting_table_S1.tab file) and save it to the .ephemeris_data directory created in your home directory. You can then run the philpott_list.py file to prepare the list for comparison with the Hollman encounters later.

Either of the directories, "encounter" or "crossings", can now be ran depending on the desire of the user. A summary of the files in each directory is below.

The general structure is that each directory is ordered by complexity. That is, the directory labelled one is conceptually the simplest analysis, while higher numbered directories are more complex.

## Crossings Directory

The crossings directory contains the files of the analysis of the Hollman et al. crossing list 2025. They can be ran in any particular order, and their titles should give a decent idea of their content.

### Generic Analysis & Visualisation

This directory looks at the number of crossings per orbit, and breaks them up by type (Magnetopause in/out & bow shock in/out). It also calculates the time between each crossing and each crossing of the same type. This temporal data contains information about highly variable events when the time between crossings is small. It also contains the distribution of times between encounters, making it a useful comparison.

### Extrema & Outlier Analysis

This takes the data from directory one and looks at the outliers of the data. It looks at the orbits with the highest number of crossings, as well as all the orbits with an odd number of crossings. It also plots the 5 crossings in the Hollman list that are labelled as unphysical.

### Crossings Close to Apoapsis & Periapsis
By looking at crossings near apoapsis (for the bow shock) and periapsis (for the magnetopause), we can see when these boundaries are highly compressed.

### Magnetosheath traversals
Using the Hollman crossings list, we take the last bow shock in and the subsequent magnetopause in to define an inbound magnetosheath traversal. Similarly, we define an outbound magnetosheath traversal as the time from the last magnetopause out to the first bow shock out. We first calculate the length of each traversal and calculate an estimate of how long MESSENGER was in the magnetosheath. The interest lies in the shortest traversals where we expect either the magnetopause or the bow shock to be compressed. We then look at the number of traversals and the length of the shortest traversals, spatially, normalised by MESSENGER's residence. The grazing angle, which is defined with respect to the normal of the average boundary position from Winslow et al., is calculated for each respective crossing. When the grazing angle is 0, it is perpendicular to the boundary, and 90 when it is parallel.

## Encounters Directory

### Generic Analysis & Visualisation
Repeats analysis above for encounters. An encounter is defined using the trajectory direction from the Hollman list. Intuitively, an encounter is a collection of crossings of the same type, from first time into the boundary, to the last time, i.e first bow shock in to last bow shock in, without being interrupted by a magnetopause crossing. This is defined in code by grouping crossings of the same type (bow shock - bow shock & magnetopause-magnetopause) that are also in the same _direction_. 

The direction is defined in the Hollman list by whether or not the velocity vector of MESSENGER dotted with the normal vector to the average boundary is positive or negative. With this definition of encounters, we group the crossings into a new list by encounters, and we can analyse the number of encounters per orbit, number of crossings per encounter, and time between encounters.

### Extrema & Outlier Analysis
This investigates the encounters with an even number of crossings. As a complete encounter should ultimately be going either in or out of a boundary, it must be opened and terminated by a crossing of the same label. Therefore, the number of crossings must be odd.

We also check the number of encounters that overlap the orbits. By summing the number of encounters found per orbit, we can see the encounters that were _overcounted_ due to the encounter lasting into the next orbit.

### Spatial Histograms
This creates 2D histograms superimposed on top of the spatial slices of MESSENGER's orbit. We plot residence, encounters normalised by residence, and number of crossings per encounter.
