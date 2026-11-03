from pathlib import Path

import pyvista as pv


def add_trajectories(
    plotter,
    vtk_file,
    color="red",
    line_width=3,
    opacity=1.0,
    verbose=True,
):
    """
    Add trajectory lines if the file exists.

    Returns
    -------
    actor or None
    """

    vtk_file = Path(vtk_file)

    if not vtk_file.exists():

        if verbose:
            print(
                f"[INFO] Trajectory file "
                f"'{vtk_file}' not found. Skipping."
            )

        return None

    try:

        trajectories = pv.read(str(vtk_file))

    except Exception as e:

        if verbose:
            print(
                f"[WARNING] Failed to read "
                f"'{vtk_file}': {e}"
            )

        return None

    if trajectories.n_points == 0:

        if verbose:
            print(
                f"[INFO] Empty trajectory file "
                f"'{vtk_file}'. Skipping."
            )

        return None

    return plotter.add_mesh(
        trajectories,
        color=color,
        line_width=line_width,
        opacity=opacity,
        render_lines_as_tubes=True,
    )