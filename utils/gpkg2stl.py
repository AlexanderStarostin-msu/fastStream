#!/usr/bin/env python3
"""
Convert OSM GeoPackage building footprints into watertight STL solids for OpenFOAM.

Height priority:
    1. OSM ``height`` tag
    2. OSM ``building:levels`` tag multiplied by 3.2 m
    3. A 10 m default

USAGE:
    python gpkg2stl.py INPUT.gpkg OUTPUT.stl

REQUIREMENTS:
    Python 3.9+ with geopandas, shapely 2.1+, trimesh, numpy, and pyvista.
    Input must contain a ``multipolygons`` layer with a ``building`` column
    and a defined coordinate reference system. Buildings are reprojected to
    EPSG:25832, then translated so the footprint bounds are centered on the
    local XY origin before geometry processing and extrusion. The ground plane
    remains at z=0. Polygon rings are topology-preserving
    simplified by 30 cm to remove redundant near-collinear vertices. Each
    extruded building is checked for watertightness and consistent
    winding. The combined mesh and reloaded STL are checked for finite
    coordinates, positive-area triangles, and a minimum normalized triangle
    quality of 0.01. Constrained Delaunay triangulation avoids sliver cap
    facets without altering footprint geometry. Overlapping footprints are
    assigned to the tallest building, and footprints are inset by 25 cm to
    prevent shared walls between separate solids. Local coordinates avoid
    precision loss in downstream mesh generation and CFD calculations.
"""

import argparse
import math
import re
from pathlib import Path

import geopandas as gpd
import numpy as np
import pyvista as pv
import trimesh
from shapely import constrained_delaunay_triangles, make_valid
from shapely.geometry import GeometryCollection, MultiPolygon, Polygon

DEFAULT_HEIGHT_M = 10.0
LEVEL_HEIGHT_M = 3.2
GEOMETRY_SIMPLIFY_TOLERANCE_M = 0.3
MIN_TRIANGLE_QUALITY = 0.01
BUILDING_SEPARATION_M = 0.25


def parse_height(tags):
    """Parse an OSM height tag in metres."""
    if tags is None:
        return None

    match = re.search(r'"height"=>"([^"]+)"', str(tags))
    if not match:
        return None

    value = match.group(1).strip().lower().removesuffix("m").strip()
    try:
        height = float(value)
    except ValueError:
        return None
    return height if math.isfinite(height) and height > 0 else None


def parse_levels(tags):
    """Parse an OSM building:levels tag."""
    if tags is None:
        return None

    match = re.search(r'"building:levels"=>"([^"]+)"', str(tags))
    if not match:
        return None

    try:
        levels = float(match.group(1))
    except ValueError:
        return None
    return levels if math.isfinite(levels) and levels > 0 else None


def determine_height(row):
    """Return a positive building height and the tag used to determine it."""
    tags = row.get("other_tags")
    height = parse_height(tags)
    if height is not None:
        return height, "height"

    levels = parse_levels(tags)
    if levels is not None:
        return levels * LEVEL_HEIGHT_M, "levels"

    return DEFAULT_HEIGHT_M, "default"


def extract_buildings(gpkg_file):
    buildings = gpd.read_file(gpkg_file, layer="multipolygons")
    if "building" not in buildings.columns:
        raise ValueError("GeoPackage multipolygons layer has no 'building' column")
    if buildings.crs is None:
        raise ValueError("GeoPackage multipolygons layer has no CRS")

    buildings = buildings[buildings["building"].notna()].copy()
    if buildings.empty:
        raise ValueError("No buildings found in the multipolygons layer")

    print(f"Found {len(buildings)} building footprints")
    return buildings


def project_to_metric(buildings):
    print("Reprojecting to EPSG:25832 ...")
    return buildings.to_crs(25832)


def center_buildings(buildings):
    """Translate projected footprints to a local frame centered in XY."""
    bounds = buildings.total_bounds
    if not np.isfinite(bounds).all():
        raise ValueError("Projected building footprints have non-finite bounds")

    xmin, ymin, xmax, ymax = bounds
    center_x = xmin + (xmax - xmin) / 2
    center_y = ymin + (ymax - ymin) / 2
    print(
        "Translating footprints to local XY origin "
        f"(projected center: {center_x:.3f}, {center_y:.3f} m)"
    )

    buildings = buildings.copy()
    buildings.geometry = buildings.geometry.translate(
        xoff=-center_x,
        yoff=-center_y,
    )
    return buildings.set_crs(None, allow_override=True)


def polygon_parts(geometry):
    """Yield polygon parts from polygonal geometries and geometry collections."""
    if isinstance(geometry, Polygon):
        yield geometry
    elif isinstance(geometry, (MultiPolygon, GeometryCollection)):
        for part in geometry.geoms:
            yield from polygon_parts(part)


def validate_mesh(mesh, description, require_watertight=True):
    if len(mesh.faces) == 0:
        raise ValueError(f"{description} has no faces")
    if not np.isfinite(mesh.vertices).all():
        raise ValueError(f"{description} contains non-finite vertex coordinates")

    zero_area_faces = np.flatnonzero(~np.isfinite(mesh.area_faces) | (mesh.area_faces <= 0))
    if len(zero_area_faces):
        sample = ", ".join(map(str, zero_area_faces[:10]))
        raise ValueError(
            f"{description} contains {len(zero_area_faces)} zero/invalid-area "
            f"triangles (first face indices: {sample})"
        )

    triangles = mesh.vertices[mesh.faces]
    edge_lengths = np.stack(
        (
            np.linalg.norm(triangles[:, 1] - triangles[:, 0], axis=1),
            np.linalg.norm(triangles[:, 2] - triangles[:, 1], axis=1),
            np.linalg.norm(triangles[:, 0] - triangles[:, 2], axis=1),
        ),
        axis=1,
    )
    qualities = 4 * math.sqrt(3) * mesh.area_faces / np.square(edge_lengths).sum(axis=1)
    poor_quality_faces = np.flatnonzero(qualities < MIN_TRIANGLE_QUALITY)
    if len(poor_quality_faces):
        sample = ", ".join(map(str, poor_quality_faces[:10]))
        raise ValueError(
            f"{description} contains {len(poor_quality_faces)} sliver triangles "
            f"with normalized quality below {MIN_TRIANGLE_QUALITY:g} "
            f"(minimum {qualities.min():g}; first face indices: {sample})"
        )

    if require_watertight and not mesh.is_watertight:
        raise ValueError(f"{description} is not watertight")
    if not mesh.is_winding_consistent:
        raise ValueError(f"{description} has inconsistent face winding")


def extrude_polygon(polygon, height):
    """Extrude a polygon with constrained caps and explicit boundary walls."""
    triangles = list(constrained_delaunay_triangles(polygon).geoms)
    if not triangles:
        raise ValueError("polygon triangulation produced no triangles")

    vertices = []
    vertex_indices = {}
    faces = []
    edge_directions = {}

    def vertex_index(coordinate, z):
        point = (float(coordinate[0]), float(coordinate[1]), z)
        if point not in vertex_indices:
            vertex_indices[point] = len(vertices)
            vertices.append(point)
        return vertex_indices[point]

    for triangle in triangles:
        coordinates = list(triangle.exterior.coords)[:3]
        a, b, c = coordinates
        signed_area = (
            (b[0] - a[0]) * (c[1] - a[1])
            - (b[1] - a[1]) * (c[0] - a[0])
        )
        if signed_area < 0:
            b, c = c, b

        top = [vertex_index(point, height) for point in (a, b, c)]
        bottom = [vertex_index(point, 0.0) for point in (a, b, c)]
        faces.extend((top, list(reversed(bottom))))
        for start, end in ((a, b), (b, c), (c, a)):
            start = (float(start[0]), float(start[1]))
            end = (float(end[0]), float(end[1]))
            key = tuple(sorted((start, end)))
            edge_directions.setdefault(key, []).append((start, end))

    triangulated_area = sum(triangle.area for triangle in triangles)
    if not math.isclose(
        triangulated_area,
        polygon.area,
        rel_tol=1e-10,
        abs_tol=1e-8,
    ):
        raise ValueError(
            "constrained triangulation does not cover the polygon area "
            f"({triangulated_area} != {polygon.area})"
        )

    for edges in edge_directions.values():
        if len(edges) == 2:
            if edges[0] == edges[1]:
                raise ValueError("cap triangulation has inconsistent edge winding")
            continue
        if len(edges) != 1:
            raise ValueError("cap triangulation has a non-manifold boundary")

        start, end = edges[0]
        bottom_a = vertex_index(start, 0.0)
        bottom_b = vertex_index(end, 0.0)
        top_a = vertex_index(start, height)
        top_b = vertex_index(end, height)
        faces.extend(
            (
                (bottom_a, bottom_b, top_b),
                (bottom_a, top_b, top_a),
            )
        )

    mesh = trimesh.Trimesh(
        vertices=np.asarray(vertices, dtype=float),
        faces=np.asarray(faces, dtype=np.int64),
        process=True,
    )
    return mesh


def build_mesh(buildings):
    meshes = []
    height_sources = {"height": 0, "levels": 0, "default": 0}
    repaired_count = 0
    records = []

    for feature_id, row in buildings.iterrows():
        height, source = determine_height(row)
        height_sources[source] += 1

        geometry = row.geometry
        if geometry is None or geometry.is_empty:
            raise ValueError(f"Building feature {feature_id} has empty geometry")
        if not geometry.is_valid:
            geometry = make_valid(geometry)
            repaired_count += 1

        parts = list(polygon_parts(geometry))
        if not parts:
            raise ValueError(
                f"Building feature {feature_id} has no polygonal geometry after repair"
            )

        records.extend((feature_id, height, polygon) for polygon in parts)

    occupied_footprints = None
    for feature_id, height, polygon in sorted(
        records,
        key=lambda record: (-record[1], str(record[0])),
    ):
        if occupied_footprints is not None:
            polygon = polygon.difference(occupied_footprints)
        occupied_footprints = (
            polygon if occupied_footprints is None
            else occupied_footprints.union(polygon)
        )
        polygon = polygon.buffer(-BUILDING_SEPARATION_M)
        parts = list(polygon_parts(polygon))
        if not parts:
            continue

        for part_index, polygon_part in enumerate(parts):
            polygon_part = polygon_part.simplify(
                GEOMETRY_SIMPLIFY_TOLERANCE_M,
                preserve_topology=True,
            )
            if (
                polygon_part.is_empty
                or not polygon_part.is_valid
                or polygon_part.area <= 0
            ):
                raise ValueError(
                    f"Building feature {feature_id}, polygon {part_index} is "
                    "empty, invalid, or has zero area after cleanup"
                )

            try:
                mesh = extrude_polygon(polygon_part, height)
                validate_mesh(
                    mesh,
                    f"Building feature {feature_id}, polygon {part_index}",
                )
            except Exception as exc:
                raise RuntimeError(
                    f"Failed to create a valid solid for building feature "
                    f"{feature_id}, polygon {part_index}: {exc}"
                ) from exc

            meshes.append(mesh)

    city = trimesh.util.concatenate(meshes)
    validate_mesh(city, "Combined building mesh")

    print("\nHeight statistics:")
    print(f"  from 'height'          : {height_sources['height']}")
    print(f"  from 'building:levels' : {height_sources['levels']}")
    print(f"  default {DEFAULT_HEIGHT_M:.1f} m      : {height_sources['default']}")
    print(f"  repaired geometries    : {repaired_count}")
    print(f"  footprint gap          : {BUILDING_SEPARATION_M:.3f} m")
    print(f"  polygon solids         : {len(meshes)}")
    return city


def export_stl(mesh, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    mesh.export(output_path, file_type="stl_ascii")

    exported_mesh = trimesh.load(output_path, force="mesh", process=False)
    validate_mesh(
        exported_mesh,
        f"Exported STL {output_path}",
        require_watertight=False,
    )
    return exported_mesh


def save_preview(mesh, output_path):
    plotter = pv.Plotter(off_screen=True)
    plotter.add_mesh(mesh, color="lightgray", show_edges=False)
    plotter.add_axes()
    plotter.camera_position = "iso"
    plotter.show(screenshot=str(output_path))
    print(f"Saved {output_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("USAGE:")[0].strip())
    parser.add_argument("input_gpkg", type=Path, help="input OSM GeoPackage")
    parser.add_argument("output_stl", type=Path, help="output STL path")
    args = parser.parse_args()

    buildings = extract_buildings(args.input_gpkg)
    buildings = project_to_metric(buildings)
    buildings = center_buildings(buildings)
    city = build_mesh(buildings)

    print(f"\nVertices: {len(city.vertices)}")
    print(f"Faces   : {len(city.faces)}")

    exported_mesh = export_stl(city, args.output_stl)
    print(
        f"\nWritten and validated: {args.output_stl} "
        f"({len(exported_mesh.faces)} triangles)"
    )

    preview_path = args.output_stl.with_name(
        f"{args.output_stl.stem}IsoView.png"
    )
    save_preview(exported_mesh, preview_path)


if __name__ == "__main__":
    main()
