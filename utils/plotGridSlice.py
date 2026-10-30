#!/usr/bin/env python3

"""Render the latest OpenFOAM mesh slice at the midpoint of its flight band.

USAGE:
    python plotGridSlice.py [CASE] [-o OUTPUT.png]

REQUIREMENTS:
    Python 3 with PyVista and Shapely, plus an OpenFOAM case containing
    system/fsStreamLineDict, a mesh under constant/polyMesh, fresh VTK output,
    and a building surface under constant/geometry. STL surfaces are
    preferred; if none exist, an OBJ surface compressed as .obj.gz is
    supported.
"""

import argparse
import gzip
import re
import shutil
import tempfile
from pathlib import Path
from typing import Optional

import pyvista as pv
from shapely import LineString, covers, points
from shapely.ops import polygonize, unary_union


def read_scalar(
    dictionary: str,
    key: str,
    default: Optional[float] = None,
) -> float:
    match = re.search(
        rf"^\s*{re.escape(key)}\s+"
        r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)\s*;",
        dictionary,
        re.MULTILINE,
    )
    if match is None:
        if default is not None:
            return default
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
        raise ValueError("the building surface has no contours at the slice height")

    footprints = list(polygonize(unary_union(line_segments)))
    if not footprints:
        raise ValueError("could not form closed building footprints from STL")
    return unary_union(footprints)


def load_building_mesh(geometry_dir: Path):
    stl_files = sorted(path for path in geometry_dir.glob("*.stl") if path.is_file())
    if stl_files:
        return pv.read(stl_files[0]), stl_files[0]

    obj_files = sorted(path for path in geometry_dir.glob("*.obj.gz") if path.is_file())
    if not obj_files:
        raise FileNotFoundError(
            f"no STL or compressed OBJ building surface found in {geometry_dir}"
        )

    compressed_obj = obj_files[0]
    with tempfile.TemporaryDirectory(prefix="plotGridSlice-") as temporary_dir:
        obj_path = Path(temporary_dir) / compressed_obj.name[:-3]
        with gzip.open(compressed_obj, "rb") as source, obj_path.open("wb") as target:
            shutil.copyfileobj(source, target)
        return pv.read(obj_path), compressed_obj


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
            "red using minObstacleClearance."
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
        clearance = read_scalar(dictionary, "minObstacleClearance", default=0)
    except (OSError, ValueError) as error:
        parser.error(str(error))

    if min_height > max_height:
        parser.error("minHeight must be less than or equal to maxHeight")
    if clearance < 0:
        parser.error("minObstacleClearance must be non-negative")
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
        key=lambda path: path.stat().st_mtime_ns,
        default=None,
    )
    if latest_vtk is None:
        parser.error(f"no VTK timestep files found for case {case_name!r}")

    mesh_files = [
        case_dir / "constant" / "polyMesh" / name
        for name in ("points", "faces")
    ]
    missing_mesh_files = [path for path in mesh_files if not path.is_file()]
    if missing_mesh_files:
        parser.error(
            "OpenFOAM mesh files not found: "
            + ", ".join(str(path) for path in missing_mesh_files)
        )
    newest_mesh_mtime = max(path.stat().st_mtime_ns for path in mesh_files)
    if latest_vtk.stat().st_mtime_ns < newest_mesh_mtime:
        parser.error(
            f"VTK mesh is stale: {latest_vtk} predates constant/polyMesh. "
            "Regenerate it with 'foamToVTK -constant -noPointValues', "
            "then rerun plotGridSlice.py."
        )

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

    geometry_dir = case_dir / "constant" / "geometry"
    try:
        building_mesh, building_surface = load_building_mesh(geometry_dir)
    except OSError as error:
        parser.error(str(error))
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
    print(f"Building surface: {building_surface}")
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
