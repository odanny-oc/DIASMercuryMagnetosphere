from astropy.table import QTable, vstack
from astropy.time import Time

import numpy as np
import datetime as dt
import os
from hermpymod.functions.ephemeris_downsampler import parse_crossing_list, parse_periapsis_data


home_dir = os.getenv("HOME")
from hermpymod.paths import DATA_DIR
data_dir = DATA_DIR
os.makedirs(data_dir, exist_ok = True)

crossing_data = parse_crossing_list()

sig_figs = 4

def round_datetime_to_second(t):
    return t.replace(microsecond=0) + dt.timedelta(seconds=round(t.microsecond / 1e6))



def initialise_encounters(crossing_data=crossing_data):
    """
    From list of crossings with direction, construct encounters list by matching crossings of the same type and 'trajectory direction', crossing by crossing. For each crossing in the list it checks to see if the "Label" (BS or MP) is the same as the next, it then checks to see if the "Trajectory Direction" is the same and that they are less than a full orbit apart (12 hours). If they are, the crossing is added to the list and the proceeding crossing and the next proceeding crossing are checked until the conditions fail. When the conditions fail the resulting list is saved and then cleared and the process repeats until we have gone through all crossings in the list.
    

    crosssing_data: Table or QTable with columns

    "UTC": Strings in the form YYYY-MM-DD HH:MM:SS
    "Label": BS_IN, BS_OUT, MP_IN, MP_OUT
    "Trajectory Direction": Inbound or Outbound
    
    Returns list of tables where every list entry is an encounter and every table entry is a crossing in that encounter.
    """
    crossing_time = Time(crossing_data["UTC"]).to_datetime()
    crossing_label = crossing_data["Label"]
    crossing_direction = crossing_data["Trajectory Direction"]

    # encounters_data is the returned list while, cross_number is the tempory list to save the encounter currently being built

    encounters_data = []
    cross_number = []

    for i in range(len(crossing_time)):
        crossing_type = crossing_label[i][:2]

        try:
            time_delta = crossing_time[i+1] - crossing_time[i]

            # Condition stipulated above. Crossings 13 hours apart are more than a full orbit away so these are discarded. This can happen with extended gaps in the data.
            if crossing_type in crossing_label[i+1] and crossing_direction[i] == crossing_direction[i+1] and time_delta<=dt.timedelta(hours=13):
                cross_number.append(crossing_data[i])

            # When condition fails, save current item and add list to encounters_data, then reset the cross_number list
            else:
                cross_number.append(crossing_data[i])
                rows = vstack(cross_number)

                encounters_data.append(rows)
                cross_number = []

        # IndexError to handle the final item in the crossing list
        except IndexError:
            cross_number.append(crossing_data[i])
            rows = vstack(cross_number)

            encounters_data.append(rows)
            cross_number = []

    return encounters_data


def encounter_finder(crossing_data=crossing_data, verbose=False, include_outliers=False):
    """
    Post-processes the raw encounter list from initialise_encounters to fix
    encounters with an even number of crossings.

    An encounter should normally start and end on the same crossing type
    (e.g. BS_IN ... BS_IN), giving an odd number of crossings. An even-length
    encounter means a boundary crossing got split across two "encounters" by
    initialise_encounters (usually because a gap or direction change broke
    the sequence early). This function detects adjacent even encounters and
    shifts one crossing from the second into the first, making both odd.

    crossing_data: Table or QTable with columns
        "UTC": Strings in the form YYYY-MM-DD HH:MM:SS
        "Label": BS_IN, BS_OUT, MP_IN, MP_OUT
        "Trajectory Direction": Inbound or Outbound

    verbose: if True, also return the list of indices that were altered
        during post-processing (useful for debugging/inspection, automatically sets include_outliers to True as to preserve index structure).
    include_outliers: if True, keep encounters that are still even-length
        after post-processing (e.g. ones with no fixable neighbour).
        If False (default), these are dropped from the returned list.

    Returns:
        encounters_post_proccessed: list of tables, one per encounter, each
            containing the individual crossings that make up that encounter.
        altered_indices (only if verbose=True): list of indices into the
            original encounter_crossings list that were shifted/merged.
    """

    if verbose:
        include_outliers = True

    # Organises crossings by pure direction and type
    encounter_crossings = initialise_encounters(crossing_data)

    encounters_post_proccessed = []
    altered_indices = []

    # Post proccessing list, main goal is to shift encounters with an even number of crossings
    for idx, encounter in enumerate(encounter_crossings):
        try:
            # Checks if encounter doesn't start and end with same crossing (means even number of crossings)
            if encounter["Label"][0] != encounter["Label"][-1]:
                # Then checks if next encounter is also even, all even encounters are BS encounters, meaning this has no risk of mixing types. If the next encounter is even, it shifts the first crossing of the next encounter back, making them both odd.
                if encounter_crossings[idx + 1]["Label"][0] != encounter_crossings[idx + 1]["Label"][-1]:
                    # Move the first crossing of the next encounter onto the end of this one
                    encounter_shifted = vstack([encounter_crossings[idx], encounter_crossings[idx + 1][0]])
                    # Drop that crossing from the next encounter
                    encounter_removed = encounter_crossings[idx + 1][1:]

                    encounters_post_proccessed.append(encounter_shifted)
                    encounters_post_proccessed.append(encounter_removed)

                    altered_indices.append(idx)

                # Checks to see if last encounter was even, meaning this encounter had already been accounted for and can be skipped.
                elif encounter_crossings[idx - 1]["Label"][0] != encounter_crossings[idx - 1]["Label"][-1]:
                    continue

                else:
                    # Even encounter with no fixable neighbour; kept as-is (will be filtered out below unless include_outliers)
                    encounters_post_proccessed.append(encounter)

            else:
                # Odd (well-formed) encounter, keep unchanged
                encounters_post_proccessed.append(encounter)

        # IndexError to handle the final item in the list, which has no "next" encounter to check against
        except IndexError:
            encounters_post_proccessed.append(encounter)
            continue

    if not include_outliers:
        # Remove any remaining even-length encounters that couldn't be fixed (14)
        encounters_post_proccessed = [i for i in encounters_post_proccessed if len(i) % 2 == 1]

    if verbose:
        return encounters_post_proccessed, altered_indices
    else:
        return encounters_post_proccessed


def parse_encounters_list(force_rebuild=False, verbose=False):
    """
    Creates/loads list of encounter times and durations with orbit number.

    Checks disk for a cached CSV of encounter summary data (one row per
    encounter, with start/end time, label, duration, and orbit number). If
    not found (or force_rebuild=True), rebuilds it from encounter_finder and
    writes it to disk for future calls.

    force_rebuild: if True, deletes any existing cached CSV first, forcing
        the encounter list to be rebuilt from scratch.
    verbose: if True, also runs encounter_finder in verbose mode (forcing a
        rebuild) and returns the list of altered indices alongside the data.

    Returns:
        encounters_data: QTable with columns
            "Time Start": UTC start time of the encounter
            "Time End": UTC end time of the encounter
            "Label": encounter type + direction, e.g. "BSI", "MPO"
            "Encounter Duration": duration in hours
            "Orbit Number": orbit number the encounter falls within
        altered_indices (only if verbose=True): indices altered during
            post-processing in encounter_finder.
    """
    # verbose mode needs indices from a fresh run of encounter_finder, so force a rebuild
    if verbose:
        force_rebuild = True

    encounters_list_dir = os.path.join(data_dir, 'hollman_encounters_list_2026.csv')

    if force_rebuild:
        try:
            os.remove(encounters_list_dir)
        except FileNotFoundError:
            pass

    # Get encounters if not already downloaded
    try:
        encounters_data = QTable.read(encounters_list_dir)
    except FileNotFoundError:
        print("File not found. Building encounter list")

        if verbose:
            encounters_list, altered_indices = encounter_finder(verbose=True)
        else:
            encounters_list = encounter_finder()

        # Start/end UTC times for each encounter
        time_start = [i["UTC"][0] for i in encounters_list]
        time_end = [i["UTC"][-1] for i in encounters_list]

        # Match each encounter's start time to the nearest preceding periapsis to get its orbit number
        orbits = parse_periapsis_data()
        orbit_indices = np.searchsorted(Time(orbits["UTC"]), Time(time_start))
        orbit_number = orbits["Orbit Number"][orbit_indices]

        # Build a short label per encounter, e.g. "BS" + "I" -> "BSI"
        label = []
        for encounter in encounters_list:
            encounter_type = encounter["Label"][0][:2]
            direction = encounter["Trajectory Direction"][0][0].upper()
            label.append(encounter_type + direction)

        # Duration of each encounter in hours
        encounter_duration = np.array([(Time(time_end[i]).to_datetime() - Time(time_start[i]).to_datetime()).total_seconds()/3600 for i in range(len(time_start))])

        time_start = [round_datetime_to_second(t) for t in Time(time_start).to_datetime()]
        time_end = [round_datetime_to_second(t) for t in Time(time_end).to_datetime()]


        encounters_data = QTable({
                "Time Start": time_start,
                "Time End": time_end,
                "Label": label,
                "Encounter Duration": np.round(encounter_duration, sig_figs,),
                "Orbit Number": orbit_number
                })
        encounters_data.write(encounters_list_dir)
        print("Saved encounter list to ", encounters_list_dir)

    if verbose:
        return encounters_data, altered_indices
    else:
        return encounters_data
