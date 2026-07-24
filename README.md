# DIAS Magnetosphere Group Computing Research Internship

This repository details the progress made on the analysis of the magnetic boundary crossings made by the MESSENGER spacecraft from 2011 to 2015, for the DIAS Computing Research Internship.

To start, run the orbit_times_analysis file. This reads the "peaks_data.csv", which is a list of times when MESSENGER was at periapsis. This file verifies that all the orbits are the correct times as well as investigates the outliers of the orbits, namely the eleven ~9.5hr transition orbits that occur due to MESSENGER's two-stage transition from 12hr to 8hr.

Next, one runs the sun_data.py and philpott_list.py files to download the alternative encounter data sets for comparison later. To run philpott_list.py, one must download the Philpott et al. list of encounters https://borealisdata.ca/dataset.xhtml?persistentId=doi:10.5683/SP2/1U6FEO (the supporting_table_S1.tab file) and save it to the .ephemeris_data directory created in your home directory. You can then run the philpott_list.py file to prepare the list for comparison with the Hollman encounters later.

Either of the directories, "encounter" or "crossings", can now be ran depending on the desire of the user. A summary of the files in each directory is below.

## crossings Directory

The crossings directory contains the files of the analysis of the Hollman et al. crossing list 2025. They can be ran in any particular order, and their titles should give a decent idea of their content.

### total_per_orbit

This plots a histogram of the total number of crossings per orbit.

### per_orbit_by_type

This further breaks down the above plot into magnetopause and bow shock crossings per orbit.

### dt_between_crossings

This creates a histogram of the time between each crossing, the time between each bow shock crossing, and the time between each magnetopause crossing.


### extrema_analysis

The rest of the directory is dedicated to analysing the outliers of the crossing histogram above. The extrema are the 42 and 54 BS crossings that happen in a single orbit.


### odd_crossing

This looks at all the odd-numbered crossings, which we expect to be unphysical. It plots the MAG data and shows that these occur due to errors or missed crossings in the Hollman list or data gaps in the MESSENGER MAG data.


### unphysical

This looks at the 5 unphysical labelled crossings in the Hollman list.

## encounter Directory

### encounters_histogram
This plots the histograms from the crossings directory, but for the encounters, as well as comparing the number of crossings per encounter.

### dt_between_encounters
 Again, this is similar to dt_between_crossings but for the encounters. To run this its crucial to run the Sun and Philpott data files from earlier, as it compares the results to those lists.

### even_encounters

This looks at a few examples of encounters with an even number of crossings.
