#!/usr/bin/env python3

"""Render the latest OpenFOAM mesh slice at the midpoint of its flight band."""

import argparse
import re
from pathlib import Path

import pyvista as pv
from shapely import LineString, covers, points
from shapely.ops import polygonize, unary_union


def read_scalar(dictionary: str, key: str) -> float:
    match = re.search(
        rf"^\s*{re.escape(key)}\s+"
        r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)\s*;",
        dictionary,
        re.MULTILINE,
    )
    if match is None:
        raise ValueError(f"{key} is missing from system/fsStreamLineDict")
    return float(match.group(1))


def building_footprints(section: pv.PolyData):
    line_segments = []
    connectivity = section.lines
    offset = 0
    while offset < len(connectivity):
        point_count = int(connectivity[offset])
        point_ids = connectivity[offset + 1 : offset + point_count + 1]
        offset += point_count + 1
        if point_count >= 2:
            line_segments.append(LineString(section.points[point_ids, :2]))

    if not line_segments:
        raise ValueError("the building STL has no contours at the slice height")

    footprints = list(polygonize(unary_union(line_segments)))
    if not footprints:
        raise ValueError("could not form closed building footprints from STL")
    return unary_union(footprints)


def add_outline(plotter: pv.Plotter, geometry, z: float, color: str) -> None:
    boundaries = [geometry] if geometry.geom_type == "Polygon" else geometry.geoms
    for boundary in boundaries:
        for ring in [boundary.exterior, *boundary.interiors]:
            coordinates = [(x, y, z + 0.02) for x, y in ring.coords]
            outline = pv.lines_from_points(coordinates, close=False)
            plotter.add_mesh(
                outline,
                color=color,
                line_width=2,
                render_lines_as_tubes=True,
            )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Save a top-down PNG of the latest VTK mesh cut at the midpoint "
            "of minHeight and maxHeight. Mesh cells near buildings are marked "
            "red using visualizationClearance."
        )
    )
    parser.add_argument(
        "case",
        nargs="?",
        type=Path,
        default=Path.cwd(),
        help="OpenFOAM case directory (defaults to the current directory)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="PNG path (defaults to <case>/<case>_grid_slice_z_<height>m.png)",
    )
    args = parser.parse_args()

    case_dir = args.case.expanduser().resolve()
    dictionary_path = case_dir / "system" / "fsStreamLineDict"
    if not dictionary_path.is_file():
        parser.error(f"configuration dictionary not found: {dictionary_path}")

    try:
        dictionary = dictionary_path.read_text()
        min_height = read_scalar(dictionary, "minHeight")
        max_height = read_scalar(dictionary, "maxHeight")
        clearance = read_scalar(dictionary, "visualizationClearance")
    except (OSError, ValueError) as error:
        parser.error(str(error))

    if min_height > max_height:
        parser.error("minHeight must be less than or equal to maxHeight")
    if clearance < 0:
        parser.error("visualizationClearance must be non-negative")
    slice_height = (min_height + max_height) / 2

    vtk_dir = case_dir / "VTK"
    if not vtk_dir.is_dir():
        parser.error(f"VTK directory not found: {vtk_dir}")

    case_name = case_dir.name
    vtk_files = [
        path
        for path in vtk_dir.glob(f"{case_name}_*.vtk")
        if path.stem.rsplit("_", 1)[-1].isdigit()
    ]
    latest_vtk = max(
        vtk_files,
        key=lambda path: int(path.stem.rsplit("_", 1)[-1]),
        default=None,
    )
    if latest_vtk is None:
        parser.error(f"no VTK timestep files found for case {case_name!r}")

    mesh = pv.read(latest_vtk)
    if not mesh.n_cells:
        parser.error(f"VTK mesh contains no cells: {latest_vtk}")

    section = mesh.slice(
        normal=(0, 0, 1),
        origin=(0, 0, slice_height),
    )
    if not section.n_cells:
        parser.error(f"mesh has no cells intersecting z={slice_height:g} m")

    output = args.output or case_dir / (
        f"{case_name}_grid_slice_z_{slice_height:g}m.png"
    )
    output = output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)

    building_stl = case_dir / "constant" / "geometry" / "LdW.stl"
    if not building_stl.is_file():
        parser.error(f"building STL not found: {building_stl}")
    building_mesh = pv.read(building_stl)
    building_section = building_mesh.slice(
        normal=(0, 0, 1),
        origin=(0, 0, slice_height),
    )
    try:
        building_geometry = building_footprints(building_section)
    except ValueError as error:
        parser.error(str(error))

    buffered_buildings = building_geometry.buffer(clearance)
    centers = section.cell_centers().points
    center_points = points(centers[:, 0], centers[:, 1])
    unsafe_cells = covers(buffered_buildings, center_points)

    pv.OFF_SCREEN = True
    plotter = pv.Plotter(off_screen=True, window_size=(2048, 1536))
    safe_mesh = section.extract_cells(~unsafe_cells)
    unsafe_mesh = section.extract_cells(unsafe_cells)
    if safe_mesh.n_cells:
        plotter.add_mesh(
            safe_mesh,
            color="#38a169",
            show_edges=True,
            edge_color="#344054",
            line_width=0.7,
            label="Motion space",
        )
    if unsafe_mesh.n_cells:
        plotter.add_mesh(
            unsafe_mesh,
            color="#e53935",
            show_edges=True,
            edge_color="#344054",
            line_width=0.7,
            label="Building / clearance",
        )

    add_outline(plotter, building_geometry, slice_height, "#6b1111")
    if clearance > 0:
        add_outline(plotter, buffered_buildings, slice_height, "#e53935")
    plotter.add_legend(bcolor="white", face="rectangle", loc="lower right")
    plotter.add_text(
        (
            f"CFD mesh slice at z = {slice_height:g} m\n"
            f"Red: building footprints + {clearance:g} m clearance"
        ),
        position="upper_left",
        font_size=12,
        color="black",
    )
    plotter.add_axes()
    plotter.view_xy()
    plotter.reset_camera()
    plotter.screenshot(str(output))
    plotter.close()

    print(f"VTK mesh: {latest_vtk}")
    print(f"Flight band: {min_height:g}-{max_height:g} m")
    print(f"Slice height: {slice_height:g} m")
    print(f"Preview building clearance: {clearance:g} m")
    print(
        f"Mesh cells: {safe_mesh.n_cells} motion-space, "
        f"{unsafe_mesh.n_cells} building/clearance"
    )
    print(f"Saved grid slice: {output}")


if __name__ == "__main__":
    main()
