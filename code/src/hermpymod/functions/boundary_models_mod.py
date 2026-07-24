import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from typing import Literal, get_args
from hermpy.utils.constants import Constants

frame_types = Literal["MSM", "MSO"]


def plot_magnetospheric_boundaries(
    ax: plt.Axes,
    plane: Literal["xy", "yz", "xz"] = "xy",
    frame: frame_types = "MSM",
    cylindrical=False,
    sub_solar_magnetopause: float = 1.45,
    alpha: float = 0.5,
    psi: float = 1.04,
    p: float = 2.75,
    initial_x: float = 0.5,
    add_legend: bool = False,
    zorder: int = 0,
    color="black",
    lw=1,
) -> None:
    """Add average magnetopause and bow shock locations based on
    Winslow et al. (2013).

    Add the plane projection of the average magnetopause and
    bow shock locations based on Winslow et al. (2013).
    These are plotted in units of Mercury radii.


    Parameters
    ----------
    ax : pyplot.Axes
        The pyplot axis to add the boundaries to.

    plane : str {`"xy"`, `"xz"`, `"yz"`}, optional
        What plane to project the boundaries to. yz is not yet
        implimented.

    add_legend : bool {`True`, `False`}, optional
        Should pyplot legend labels be added.


    Returns
    -------
    None
    """

    # Define dipole offset for conversion to MSO
    Zd = Constants.DIPOLE_OFFSET.to("Mercury Radii").value

    # Plotting magnetopause
    phi = np.linspace(0, 2 * np.pi, 1000)
    rho_mp_boundary = sub_solar_magnetopause * (2 / (1 + np.cos(phi))) ** alpha

    magnetopause_x_coords = rho_mp_boundary * np.cos(phi)
    magnetopause_y_coords = rho_mp_boundary * np.sin(phi)
    if frame == "MSM":
        magnetopause_z_coords = magnetopause_y_coords
    elif frame == "MSO":
        magnetopause_z_coords = magnetopause_y_coords + Zd
    else:
        raise ValueError(
            f"Invalid frame: {frame!r}. Must be one of {get_args(frame_types)}"
        )

    L = psi * p

    rho_bs_boundary = L / (1 + psi * np.cos(phi))

    bowshock_x_coords = initial_x + rho_bs_boundary * np.cos(phi)
    bowshock_y_coords = rho_bs_boundary * np.sin(phi)

    # Bow shock functional form creates non-physical points far sunward of Mercury.
    # These are incorrect and must be removed.
    bowshock_y_coords = bowshock_y_coords[bowshock_x_coords < 2]
    bowshock_x_coords = bowshock_x_coords[bowshock_x_coords < 2]

    if frame == "MSM":
        bowshock_z_coords = bowshock_y_coords
    elif frame == "MSO":
        bowshock_z_coords = bowshock_y_coords + Zd
    else:
        raise ValueError(
            f"Invalid frame: {frame!r}. Must be one of {get_args(frame_types)}"
        )

    plane_coordinates = {
        "xy": (
            bowshock_x_coords,
            bowshock_y_coords,
            magnetopause_x_coords,
            magnetopause_y_coords,
        ),
        "yz": (
            bowshock_y_coords,
            bowshock_z_coords,
            magnetopause_y_coords,
            magnetopause_z_coords,
        ),
        "xz": (
            bowshock_x_coords,
            bowshock_z_coords,
            magnetopause_x_coords,
            magnetopause_z_coords,
        ),
    }

    # Set coordinates for desired plane
    (
        bowshock_plot_x_coords,
        bowshock_plot_y_coords,
        magnetopause_plot_x_coords,
        magnetopause_plot_y_coords,
    ) = plane_coordinates[plane]

    bowshock_plot_label = ""
    magnetopause_plot_label = ""

    if add_legend:
        bowshock_plot_label = "Avg. Bowshock (Winslow et al. 2013)"
        magnetopause_plot_label = "Avg. Magnetopause (Winslow et al. 2013)"

    if cylindrical:
        if frame == "MSM":
            rho_mp = np.abs(magnetopause_y_coords)
            rho_bs = np.abs(bowshock_y_coords)
        elif frame=="MSO": 
            rho_mp = np.abs(magnetopause_z_coords)
            rho_bs = np.abs(bowshock_z_coords)
        else:
            raise ValueError

        ax.plot(
            magnetopause_x_coords,
            rho_mp,
            ls="--",
            lw=lw,
            color=color,
            label=magnetopause_plot_label,
            zorder=zorder,
        )
        ax.plot(
            bowshock_x_coords,
            rho_bs,
            ls="-",
            lw=lw,
            color=color,
            label=bowshock_plot_label,
            zorder=zorder,
        )

    elif plane == "yz":
        # Plot concentric circles of boundaries
        x_coords = np.arange(0,10, 0.5)
        cmap = plt.cm.plasma_r
        norm = plt.Normalize(vmin=x_coords.min(), vmax=x_coords.max())

        for x in x_coords:
            mp_radius = 2* sub_solar_magnetopause - x 
            bs_radius = L - psi*x
            if mp_radius < 0:
                pass
            else:
                ax.plot(
                    mp_radius * np.cos(phi),
                    mp_radius * np.sin(phi) + Zd,
                    alpha = 0.3,
                    ls="--",
                    lw=lw,
                    color=cmap(norm(x)),
                    label=magnetopause_plot_label,
                    zorder=zorder,
                )
            ax.plot(
                bs_radius * np.cos(phi),
                bs_radius * np.sin(phi) + Zd,
                ls="-",
                alpha = 0.3,
                lw=lw,
                color=cmap(norm(x)),
                label=bowshock_plot_label,
                zorder=zorder,
            )

    else:
        ax.plot(
            magnetopause_plot_x_coords,
            magnetopause_plot_y_coords,
            ls="--",
            lw=lw,
            color=color,
            label=magnetopause_plot_label,
            zorder=zorder,
        )
        ax.plot(
            bowshock_plot_x_coords,
            bowshock_plot_y_coords,
            ls="-",
            lw=lw,
            color=color,
            label=bowshock_plot_label,
            zorder=zorder,
        )
