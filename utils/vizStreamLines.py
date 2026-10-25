#!/usr/bin/env python3

# Requirements:
# - Read the latest timestep mesh for the selected case.
# - Seed forward streamlines from every face on inlet-like patches.
# - Color each streamline by its total time of flight and show a scalar bar.
# - Display the mesh edges in the requested coordinate-plane view, or isometric
#   when no plane is specified.
# - Save the image as <case_name>.png in the current working directory.

# Usage
#   In advance the foamToVTK should be run after the OpenFOAM case simulation
#   ./vizStreamLines.py "/workspaces/fastStream/OpenFOAM-10/tutorials/incompressible/simpleFoam/pitzDaily" xy

import argparse
from pathlib import Path

import numpy as np
import pyvista as pv

parser = argparse.ArgumentParser()
parser.add_argument(
    "case",
    nargs="?",
    help="case directory or case name (defaults to the current directory)",
)
parser.add_argument(
    "plane",
    nargs="?",
    type=str.lower,
    choices=("xz", "yz", "xy", "yx", "zy", "zx"),
    help="coordinate plane to view (defaults to isometric)",
)
args = parser.parse_args()

if args.case is None:
    case_dir = Path.cwd()
    case_name = case_dir.name
else:
    case_path = Path(args.case).expanduser()
    if case_path.is_dir():
        case_dir = case_path.resolve()
        case_name = case_dir.name
    elif case_path.is_absolute() or len(case_path.parts) > 1:
        parser.error(f"case directory does not exist: {case_path}")
    else:
        case_dir = Path.cwd()
        case_name = args.case

vtk_dir = case_dir / "VTK"
if not vtk_dir.is_dir():
    parser.error(f"VTK directory not found for case {case_name!r}: {vtk_dir}")

vtk_files = (
    path
    for path in vtk_dir.glob(f"{case_name}_*.vtk")
    if path.stem.rsplit("_", 1)[-1].isdigit()
)
latest_vtk = max(
    vtk_files,
    key=lambda path: int(path.stem.rsplit("_", 1)[-1]),
    default=None,
)
if latest_vtk is None:
    parser.error(f"no VTK timestep files found for case {case_name!r}")

mesh = pv.read(latest_vtk)

inlet_keywords = (
    "inlet",
    "front",
    "inflow",
    "enter",
    "entrance",
    "entry",
    "intake",
    "feed",
    "supply",
    "upstream",
    "injection",
    "admission",
    "source",
)
timestep = latest_vtk.stem.rsplit("_", 1)[-1]
inlet_patches = [
    patch_dir
    for patch_dir in sorted(vtk_dir.iterdir())
    if patch_dir.is_dir()
    and any(keyword in patch_dir.name.casefold() for keyword in inlet_keywords)
]
if not inlet_patches:
    parser.error(
        "no inlet-like VTK patch directories found; tried names containing "
        + ", ".join(inlet_keywords)
    )

patch_files = [
    patch_file
    for patch_dir in inlet_patches
    for patch_file in sorted(patch_dir.glob(f"*_{timestep}.vtk"))
]
patches_without_timestep = [
    patch_dir
    for patch_dir in inlet_patches
    if not any(patch_dir.glob(f"*_{timestep}.vtk"))
]
if patches_without_timestep:
    parser.error(
        f"no VTK patch file for timestep {timestep} in: "
        + ", ".join(patch_dir.name for patch_dir in patches_without_timestep)
    )

seed_points = []
for patch_file in patch_files:
    patch = pv.read(patch_file)
    seed_points.extend(patch.cell_centers().points)

if not seed_points:
    parser.error(f"no faces found on inlet-like patches at timestep {timestep}")

source = pv.PolyData(seed_points)

streamlines = mesh.streamlines_from_source(
    source,
    vectors="U",
    max_length=1000,
    integration_direction="forward",
)
if "IntegrationTime" not in streamlines.point_data:
    parser.error("streamline output does not contain integration-time values")

integration_time = streamlines.point_data["IntegrationTime"]
time_of_flight = np.array(
    [
        integration_time[streamlines.get_cell(cell_id).point_ids[-1]]
        for cell_id in range(streamlines.n_cells)
    ]
)
streamlines.cell_data["TimeOfFlight"] = time_of_flight

pv.OFF_SCREEN = True

p = pv.Plotter(off_screen=True)
p.add_mesh(
    mesh.extract_surface(algorithm="dataset_surface"),
    color="white",
    show_edges=True,
    edge_color="black",
    line_width=1,
    opacity=0.25,
)
p.add_mesh(
    streamlines,
    scalars="TimeOfFlight",
    cmap="viridis",
    line_width=3,
    scalar_bar_args={"title": "Time of flight (s)"},
)
plane_orientation = {
    "xy": ((0, 0, 1), (0, 1, 0)),
    "yx": ((0, 0, -1), (1, 0, 0)),
    "xz": ((0, -1, 0), (0, 0, 1)),
    "zx": ((0, 1, 0), (1, 0, 0)),
    "yz": ((1, 0, 0), (0, 0, 1)),
    "zy": ((-1, 0, 0), (0, 1, 0)),
}
if args.plane is None:
    p.view_isometric()
else:
    view_direction, view_up = plane_orientation[args.plane]
    p.view_vector(view_direction, viewup=view_up)
p.reset_camera()
p.screenshot(f"{case_name}.png")