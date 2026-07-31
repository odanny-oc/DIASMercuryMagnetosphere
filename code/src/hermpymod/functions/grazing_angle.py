import datetime as dt
import multiprocessing
from functools import lru_cache
from typing import Iterable, Sequence, Union, Literal, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.signal
import scipy.spatial
import spiceypy as spice
from tqdm import tqdm

from astropy.time import Time

from hermpy.utils import Constants
from hermpy.net import ClientSPICE
from hermpymod.functions.ephemeris_downsampler import parse_spice


Zd = Constants.DIPOLE_OFFSET.to("Mercury Radii")

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


def get_boundary_normal(
    position: Tuple[float, float], boundary: Literal["Bow Shock", "Magnetopause"]
) -> Tuple[float, float]:

    match boundary:
        case "Bow Shock":
            initial_x = 0.5
            psi = 1.04
            p = 2.75

            L = psi * p

            phi = np.linspace(0, 2 * np.pi, 10000)
            rho = L / (1 + psi * np.cos(phi))

            # Cylindrical coordinates (X, R)
            x_coords = initial_x + rho * np.cos(phi)
            r_coords = rho * np.sin(phi)

            boundary_positions = np.array([x_coords, r_coords]).T

        case "Magnetopause":
            sub_solar_point = 1.45
            alpha = 0.5

            phi = np.linspace(0, 2 * np.pi, 10000)
            rho = sub_solar_point * (2 / (1 + np.cos(phi))) ** alpha

            # Cylindrical coordinates (X, R)
            x_coords = rho * np.cos(phi)
            r_coords = rho * np.sin(phi)

            boundary_positions = np.array([x_coords, r_coords]).T

        case _:
            raise ValueError(
                f"Invalid function choice: {boundary}. Options are 'Bow Shock', 'Magnetopause'."
            )

    # We need to determine which point on the boundary curve is closest to the spacecraft
    # This method, utilising a k-d tree is computationally faster than iterrating through
    # the points and determining the distance.
    # O(logN) vs O(N)
    kd_tree = scipy.spatial.KDTree(boundary_positions)
    _, closest_index = kd_tree.query(position)
    closest_position = boundary_positions[closest_index]

    # We can numerically find the tangent at each point on the curve using the
    # gradient.
    dx_dphi = np.gradient(x_coords, phi)
    dr_dphi = np.gradient(r_coords, phi)

    tangent = np.array([dx_dphi[closest_index], dr_dphi[closest_index]])
    tangent = tangent / np.sqrt(np.sum(tangent**2))

    # Rotate tangent by 90 degrees to get a normal candidate.
    normal_vector = np.array([-tangent[1], tangent[0]])

    # Ensure the normal points outwards by determining the dot product between
    # the candidate normal vector and the vector to the closest position on the
    # curve.
    if np.dot(normal_vector, closest_position) < 0:
        normal_vector = -normal_vector

    return normal_vector


def get_grazing_angle(
    crossing,
    function: Literal["Bow Shock", "Magnetopause"] = "Bow Shock",
    frame: Literal["MSO", "MSM"] = "MSO",
    return_vectors: bool = False,
    verbose: bool = False,
):
    """Determine the grazing angle for a given boundary crossing

    We determine the grazing angle by comparing the velocity vector of
    MESSENGER, to the surface normal of the boundary for that crossing. The
    surface normal is determined in the following way:

    We find the closest position on the Winslow (2013) average BS and MP model
    Assuming any expansion / compression occurs parallel to the normal vector
    of the curve, the vector to the closest point on the BS / MP to MESSENGER
    is parallel with the bow shock normal at that closest point.

    These two vectors are determined in the MSM' cylindrical coordinate system
    (X MSM', sqrt( (Y MSM')^2 + (Z MSM')^2 )).


    Parameters
    ----------
    crossing :
        Crossing object as saved in Hollman et al. crossing list 20205

        Must contain columns matching:
        'UTC'
        'X MSO'
        'Y MSO'
        'Z MSO'
        'Label'

    function : str {Bow Shock, Magnetopause}
        Which boundary function to compare against

    return_vectors : bool {False, True}, optional
        Returns also the normal and velocity vectors along with the
        grazing angle.

    aberrate : bool | str {True, False, "average"}

    verbose : bool {False, True}, optional
        Prints extra information to terminal


    Returns
    -------
    grazing_angle : float
        The grazing angle in degrees for that crossing
    """

    # if isinstance(crossing, Iterable) and not isinstance(crossing, pd.Series):
    #     print("Using vectorised grazing angle calculation")
    #     return Get_Grazing_Angle_Vectorised(
    #         crossing, function, return_vectors, aberrate, verbose
    #     )

    if frame=="MSO":
        start_position = crossing["X MSO", "Y MSO", "Z MSO"]

    elif frame=="MSM":
        crossing["Z MSO"] = crossing["Z MSO"] + Zd.value
        start_position = crossing["X MSO", "Y MSO", "Z MSO"]

    next_position = parse_spice((Time(crossing["UTC"]).to_datetime() + dt.timedelta(seconds= 1)), frame=frame)

    xlab, ylab ,zlab = [i for i in next_position.keys() if frame in i]

    next_position = next_position[0]

    cylindrical_start_position = np.array(
        [start_position["X MSO"], np.sqrt(start_position["Y MSO"] ** 2 + start_position["Z MSO"] ** 2)]
    )

    cylindrical_next_position = np.array(
        [next_position[xlab].value, np.sqrt(next_position[ylab].value ** 2 + next_position[zlab].value ** 2)]
    )

    # Mercury radii/s
    cylindrical_velocity = cylindrical_next_position - cylindrical_start_position

    # normalise velocity
    cylindrical_velocity /= np.sqrt(np.sum(cylindrical_velocity**2))


    # Get the normal vector of the BS at this point
    # This is just the normalised vector between the spacecraft and the closest point,
    # as the vector between an arbitrary point and the closest point on an arbitrary
    # curve is parallel to the normal vector of that curve at that closest point.
    normal_vector = get_boundary_normal(cylindrical_start_position, boundary = function)

    # If the x component of the normal vector is negative, the vector found is
    # actually the inward pointing normal. Hence, we need to flip the vector.
    if normal_vector[0] < 0:
        normal_vector = normal_vector * -1

    grazing_angle = np.arccos(
        np.dot(normal_vector, cylindrical_velocity)
        / (np.sqrt(np.sum(normal_vector**2)) * np.sqrt(np.sum(cylindrical_velocity**2)))
    )
    
    # Radians to degrees
    grazing_angle *= 180/np.pi

    # If the grazing angle is greater than 90, then we take 180 - angle as its from the other side
    # This occurs as we don't make an assumption as to what side of the model boundary we are.
    # i.e. we could be referencing the normal, or the anti-normal.
    if grazing_angle > 90:
        # If the angle is greater than 90 degrees, we have the normal vector
        # the wrong way around. i.e. the inward pointing normal.
        grazing_angle = 180 - grazing_angle

    if verbose:
        print(f"Crossing Start Time: {crossing['Start Time']}")
        print(f"Crossing Type: {crossing['Type']}")
        print(f"Spacecraft Position: {cylindrical_start_position}")
        print(f"Normal Vector (MSM): {normal_vector}")
        print(f"Velocity Vector (MSM): {cylindrical_velocity}")
        print(f"Grazing Angle: {grazing_angle:.3f} deg.")

    if return_vectors:
        return grazing_angle, normal_vector, cylindrical_velocity

    return float(grazing_angle)

# def Get_Grazing_Angle_Vectorised(
#     crossings,
#     function: Literal["Bow Shock", "Magnetopause"] = "Bow Shock",
#     return_vectors: bool = False,
#     aberrate: bool | str = True,
# ):
#     print(f"Processing {len(crossings)} crossings")
#
#     mid_crossing_times = (
#         crossings["Start Time"] + (crossings["End Time"] - crossings["Start Time"]) / 2
#     )
#     next_times = mid_crossing_times + pd.Timedelta(seconds=1)
#
#     start_positions = (
#         np.array(
#             traj.Get_Position(
#                 "MESSENGER",
#                 mid_crossing_times,
#                 frame="MSM",
#                 aberrate=aberrate,
#             )
#         )
#         / Constants.MERCURY_RADIUS.to("km")
#     )
#
#     next_positions = (
#         np.array(
#             traj.Get_Position(
#                 "MESSENGER",
#                 next_times,
#                 frame="MSM",
#                 aberrate=aberrate,
#             )
#         )
#         / Constants.MERCURY_RADIUS.to("km")
#     )
#
#     cylindrical_start_positions = np.column_stack(
#         [
#             start_positions[:, 0],
#             np.sqrt(start_positions[:, 1] ** 2 + start_positions[:, 2] ** 2),
#         ]
#     )
#     cylindrical_next_positions = np.column_stack(
#         [
#             next_positions[:, 0],
#             np.sqrt(next_positions[:, 1] ** 2 + next_positions[:, 2] ** 2),
#         ]
#     )
#
#     cylindrical_velocities = cylindrical_next_positions - cylindrical_start_positions
#     cylindrical_velocities /= np.linalg.norm(cylindrical_velocities, axis=1)[:, None]
#
#     normal_vectors = [get_boundary_normal(position, boundary=function) for position in cylindrical_start_positions]
#
#     dot_products = np.sum(normal_vectors * cylindrical_velocities, axis=1)
#
#     grazing_angles = np.arccos(dot_products)
#     grazing_angles = np.degrees(grazing_angles)  # convert to degrees
#
#     # If the grazing angle is greater than 90, then we take 180 - angle as its from the other side
#     # This occurs as we don't make an assumption as to what side of the model boundary we are.
#     # i.e. we could be referencing the normal, or the anti-normal.
#     grazing_angles = np.where(grazing_angles > 90, 180 - grazing_angles, grazing_angles)
#
#     if return_vectors:
#         normal_vectors = np.where(grazing_angles > 90, -normal_vectors, normal_vectors)
#         return grazing_angles, normal_vectors, cylindrical_velocities
#
#     return grazing_angles
