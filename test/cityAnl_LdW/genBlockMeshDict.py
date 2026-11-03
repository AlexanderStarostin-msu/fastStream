#!/usr/bin/env python3
"""Generate a blockMeshDict sized around the city building geometry.

USAGE:
    python3 genBlockMeshDict.py

REQUIREMENTS:
    Run from the cityAnl_LdW case directory with trimesh installed and
    constant/geometry/LdW.stl already generated. Set meshDomain.horizontalOffset
    in system/snappyHexMeshDict to the desired horizontal clearance in metres.
    The script centers the STL bounds on the local XY origin, shifts
    refinementBox by the same amount, writes system/blockMeshDict with that
    offset around the building geometry, and updates locationInMesh to a point
    inside the west clearance and outside the building geometry. Run it before
    surfaceFeatures so derived feature geometry uses the shifted STL.
"""

import math
from pathlib import Path
import re

import trimesh

NUMBER_PATTERN = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"


def translate_refinement_box(dictionary, x_offset, y_offset):
    """Translate refinementBox XY bounds into the local mesh frame."""
    box_match = re.search(
        r"\brefinementBox\s*\{(?P<body>[^{}]*)\}",
        dictionary,
        re.DOTALL,
    )
    if box_match is None:
        if x_offset != 0 or y_offset != 0:
            raise ValueError(
                "system/snappyHexMeshDict must define refinementBox min/max "
                "to shift it with the building geometry"
            )
        return dictionary

    body = box_match.group("body")
    coordinate_pattern = (
        rf"\b(?P<key>min|max)\s*\(\s*({NUMBER_PATTERN})\s+"
        rf"({NUMBER_PATTERN})\s+({NUMBER_PATTERN})\s*\)\s*;"
    )
    matches = list(re.finditer(coordinate_pattern, body))
    coordinates = {match.group("key"): match for match in matches}
    if set(coordinates) != {"min", "max"}:
        raise ValueError(
            "system/snappyHexMeshDict refinementBox must define min and max"
        )

    def translate_coordinate(match):
        x, y, z = (float(match.group(index)) for index in (2, 3, 4))
        return (
            f"{match.group('key')} "
            f"({x + x_offset:.12g} {y + y_offset:.12g} {z:.12g});"
        )

    body = re.sub(coordinate_pattern, translate_coordinate, body)

    return (
        dictionary[: box_match.start("body")]
        + body
        + dictionary[box_match.end("body") :]
    )


def refinement_box_bounds(dictionary):
    """Return refinementBox XY bounds, or None when no refinementBox exists."""
    box_match = re.search(
        r"\brefinementBox\s*\{(?P<body>[^{}]*)\}",
        dictionary,
        re.DOTALL,
    )
    if box_match is None:
        return None

    coordinates = {}
    for match in re.finditer(
        rf"\b(?P<key>min|max)\s*\(\s*({NUMBER_PATTERN})\s+"
        rf"({NUMBER_PATTERN})\s+({NUMBER_PATTERN})\s*\)\s*;",
        box_match.group("body"),
    ):
        coordinates[match.group("key")] = (
            float(match.group(2)),
            float(match.group(3)),
        )
    if set(coordinates) != {"min", "max"}:
        raise ValueError(
            "system/snappyHexMeshDict refinementBox must define min and max"
        )

    xmin, ymin = coordinates["min"]
    xmax, ymax = coordinates["max"]
    bounds = (xmin, ymin, xmax, ymax)
    if not all(math.isfinite(value) for value in bounds):
        raise ValueError(
            "system/snappyHexMeshDict refinementBox has non-finite bounds"
        )
    if xmax <= xmin or ymax <= ymin:
        raise ValueError("system/snappyHexMeshDict refinementBox has empty bounds")
    return bounds


snappy_dict_path = Path("system/snappyHexMeshDict")
snappy_dict = snappy_dict_path.read_text(encoding="utf-8")
offset_match = re.search(
    r"\bmeshDomain\s*\{\s*horizontalOffset\s+"
    rf"({NUMBER_PATTERN})\s*;\s*\}}",
    snappy_dict,
)
if offset_match is None:
    raise ValueError(
        f"{snappy_dict_path} must define meshDomain.horizontalOffset"
    )
horizontal_offset = float(offset_match.group(1))
if not math.isfinite(horizontal_offset) or horizontal_offset <= 0:
    raise ValueError("meshDomain.horizontalOffset must be a positive finite value")

geometry_path = Path("constant/geometry/LdW.stl")
mesh = trimesh.load(geometry_path)

xmin, ymin, zmin = map(float, mesh.bounds[0])
xmax, ymax, zmax = map(float, mesh.bounds[1])
building_bounds = (xmin, xmax, ymin, ymax, zmin, zmax)
if not all(math.isfinite(value) for value in building_bounds):
    raise ValueError(f"{geometry_path} has non-finite bounds")
if xmax <= xmin or ymax <= ymin or zmax <= zmin:
    raise ValueError(f"{geometry_path} has empty bounds")

center_x = xmin + (xmax - xmin) / 2
center_y = ymin + (ymax - ymin) / 2
shift_x = 0 if math.isclose(center_x, 0, abs_tol=1e-9) else -center_x
shift_y = 0 if math.isclose(center_y, 0, abs_tol=1e-9) else -center_y
mesh_shifted = shift_x != 0 or shift_y != 0
if mesh_shifted:
    snappy_dict = translate_refinement_box(snappy_dict, shift_x, shift_y)
    mesh.apply_translation((shift_x, shift_y, 0))
    print(
        "Shifted building geometry to local XY origin "
        f"(translation: {shift_x:.3f}, {shift_y:.3f} m)"
    )

xmin, ymin, zmin = map(float, mesh.bounds[0])
xmax, ymax, zmax = map(float, mesh.bounds[1])
refinement_bounds = refinement_box_bounds(snappy_dict)
if refinement_bounds is not None:
    box_xmin, box_ymin, box_xmax, box_ymax = refinement_bounds
    if (
        box_xmin > xmin
        or box_ymin > ymin
        or box_xmax < xmax
        or box_ymax < ymax
    ):
        box_width = box_xmax - box_xmin
        box_height = box_ymax - box_ymin
        if box_width < xmax - xmin or box_height < ymax - ymin:
            raise ValueError(
                "system/snappyHexMeshDict refinementBox is smaller than the "
                "building geometry"
            )
        box_center_x = box_xmin + box_width / 2
        box_center_y = box_ymin + box_height / 2
        geometry_center_x = xmin + (xmax - xmin) / 2
        geometry_center_y = ymin + (ymax - ymin) / 2
        box_shift_x = geometry_center_x - box_center_x
        box_shift_y = geometry_center_y - box_center_y
        snappy_dict = translate_refinement_box(
            snappy_dict,
            box_shift_x,
            box_shift_y,
        )
        print(
            "Aligned refinementBox with local building geometry "
            f"(translation: {box_shift_x:.3f}, {box_shift_y:.3f} m)"
        )

H  = zmax - zmin

# CFD domain dimensions
top        = 6.0 * H

x0 = xmin - horizontal_offset
x1 = xmax + horizontal_offset

y0 = ymin - horizontal_offset
y1 = ymax + horizontal_offset

z0 = 0.0
z1 = top

location_match = re.search(
    rf"\blocationInMesh\s*\(\s*({NUMBER_PATTERN})\s+"
    rf"({NUMBER_PATTERN})\s+({NUMBER_PATTERN})\s*\)\s*;",
    snappy_dict,
)
if location_match is None:
    raise ValueError(f"{snappy_dict_path} must define locationInMesh")

location_in_mesh = (
    xmin - horizontal_offset / 2,
    (ymin + ymax) / 2,
    float(location_match.group(3)),
)
if not z0 < location_in_mesh[2] < z1:
    raise ValueError(
        f"locationInMesh z coordinate must be inside the generated block "
        f"domain (0, {z1:g})"
    )

snappy_dict = (
    snappy_dict[: location_match.start()]
    + "locationInMesh\n    (\n        "
    + " ".join(f"{coordinate:.12g}" for coordinate in location_in_mesh)
    + "\n    );"
    + snappy_dict[location_match.end() :]
)
if mesh_shifted:
    mesh.export(geometry_path, file_type="stl_ascii")
snappy_dict_path.write_text(snappy_dict, encoding="utf-8")

print(f"Building bounds from {geometry_path}:")
print(f"  x: {xmin:.3f} .. {xmax:.3f}")
print(f"  y: {ymin:.3f} .. {ymax:.3f}")
print(f"  z: {zmin:.3f} .. {zmax:.3f}")
print(f"Block domain bounds ({horizontal_offset:g} m horizontal offset):")
print(f"  x: {x0:.3f} .. {x1:.3f}")
print(f"  y: {y0:.3f} .. {y1:.3f}")
print(f"  z: {z0:.3f} .. {z1:.3f}")
print(
    "snappyHexMesh locationInMesh: "
    + " ".join(f"{coordinate:.3f}" for coordinate in location_in_mesh)
)

nx = 60
ny = 40
nz = 10

Path("system").mkdir(exist_ok=True)

with open("system/blockMeshDict", "w") as f:

    f.write(f"""
FoamFile
{{
    version 2.0;
    format ascii;
    class dictionary;
    object blockMeshDict;
}}

// Generated by genBlockMeshDict.py on 2026-10-25.

convertToMeters 1;

vertices
(
    ({x0} {y0} {z0})
    ({x1} {y0} {z0})
    ({x1} {y1} {z0})
    ({x0} {y1} {z0})

    ({x0} {y0} {z1})
    ({x1} {y0} {z1})
    ({x1} {y1} {z1})
    ({x0} {y1} {z1})
);

blocks
(
    hex (0 1 2 3 4 5 6 7)
    ({nx} {ny} {nz})
    simpleGrading (1 1 1)
);

edges();

boundary
(
    inlet
    {{
        type patch;
        faces ((0 4 7 3));
    }}

    outlet
    {{
        type patch;
        faces ((1 2 6 5));
    }}

    sides
    {{
        type symmetry;
        faces
        (
            (0 1 5 4)
            (3 7 6 2)
        );
    }}

    top
    {{
        type symmetry;
        faces ((4 5 6 7));
    }}

    ground
    {{
        type wall;
        faces ((0 3 2 1));
    }}
);

mergePatchPairs();
""")