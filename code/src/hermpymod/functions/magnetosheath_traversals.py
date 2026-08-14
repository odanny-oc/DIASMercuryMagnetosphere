from astropy.time import Time
from hermpymod.functions.ephemeris_downsampler import parse_crossing_list


crossing_list = parse_crossing_list()


def magnetosheath_crossings(crossing_list=crossing_list):
    """
    Function to pair the crossings just before and after the magnetosheath

    crossing_list: QTable with columns
        "UTC": Time as a string in the format YYYY-MM-DD HH:MM:SS
        "Label": Type of crossing "BS_OUT", "BS_IN", "MP_IN", "MP_OUT"


    
    Returns: ms_out and ms_in, for inbound and outbound magnetosheath traversals in the form of lists of arrays of two column objects, the crossing data for the inbound and outbound crossings.

    """

    crossing_time = Time(crossing_list["UTC"]).to_datetime()
    crossing_label = crossing_list["Label"]

    ms_out = []
    ms_in = []


    for idx, crossing in enumerate(crossing_list):
        # Checks to see the crossing type (Outbound)
        if crossing["Label"] == "MP_OUT":
            # Checks to see if its the last crossing into the magnetosheath
            if crossing_label[idx + 1] == "BS_OUT":
                ms_out.append([crossing, crossing_list[idx + 1]])
            else:
                continue
        elif crossing["Label"] == "BS_IN":
            # Same for inbound traversals
            if crossing_label[idx + 1] == "MP_IN":
                ms_in.append([crossing, crossing_list[idx + 1]])
            else:
                continue
        

    return ms_out , ms_in
