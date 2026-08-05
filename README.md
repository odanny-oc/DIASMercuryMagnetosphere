# DIAS Magnetosphere Group Computing Research Internship

This repository details the progress made on the analysis of the magnetic boundary crossings made by the MESSENGER spacecraft from 2011 to 2015, for the DIAS Computing Research Internship. The goal of this project was to analyse the extrema events of Mercury's boundaries to extract information on the compression of the boundaries despite having no upstream data for large times of the MESSENGER mission.


## Data Handling and Preparation Directory
The directories are ordered by number in order of complexity. That is, the directory labelled one is conceptually the simplest analysis, while higher numbered directories are more complex. The rule of thumb is to run the directories in numerical order, but depending on the user's purposes this is not always necessary. Directories labelled 00, are to indicate that these are not necessary to be run, but in the case of the alternate data sets, this directory needs to be run for analysis between Hollman et al., Sun et al., and Philpott et al. in the encounters directory. 


### Alternative Data Sets
The alternative data sets directory prepares the data from the Sun et al. encounters list (https://zenodo.org/records/18236915) and the Philpott et al. encounters list (https://borealisdata.ca/dataset.xhtml?persistentId=doi:10.5683/SP2/1U6FEO) for use in the later directories. It prepares the data to be in the form of 5 columns.

- Time Start
- Time End
- Label
- Encounter duration

where _Time Start_ and _Time End_ are the start and end times of an encounter with a magnetic boundary, _Label_ is the type of boundary and its direction (magnetopause/bow shock in/out, in the form of (MPO, MPI, BSO, BSI), and _Encounter Duration_ is the duration of the encounter (_Time End_ - _Time Start_) in hours.

### Orbit Finder
*Note* it is unnecessary to run any file in this directory for the analysis of crossings and encounters.

The first objective of the project is to define when an orbit starts and ends. This was done with respect to the _periapsis_. The periapsis and apoapsis are found using the scipy find_peaks function. The data set is constructed by orbit_periapsis_and_apoapsis_finder.py and is saved and downloaded from Zenodo (https://zenodo.org/records/21393912). 

The data is constructed by downsampling the SPICE kernels data (one point per minute) over the entire MESSENGER mission using the hermpy package (https://hermpy.readthedocs.io/en/stable/generated_examples/spice.html). The distance from Mercury is calculated by taking the norm of the vector in _Mercury Solar Orbiter_ (MSO) coordinates. The peak and troughs of this data are used to define apoapsis and periapsis respectively and, by extension, the orbits. The downsampled data is explicitly constructed in the downsample_spice_data python file, but is automatically constructed or called when needed.

The validity of the orbits found is verified in the orbit_times_analysis file. This reads the downloaded periapsis data file and verifies that all the orbits are the correct duration as well as investigates the outliers of the orbits, namely the eleven ~9.5hr transition orbits that occur due to MESSENGER's two-stage transition from 12hr to 8hr.

Finally, the orbit_data_plotter plots all ephemeris and MAG data for any given orbit.

Either of the directories, "encounter" or "crossings", can now be ran depending on the desire of the user. A summary of the files in each directory is below.


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


# Dependinces
All the dependincies can be installed from the requirements.txt file in the code directory. The major dependincies are

- matplotlib
- numpy
- hermpy 
- sunpy
- astropy
