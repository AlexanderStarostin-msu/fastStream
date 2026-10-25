#!/usr/bin/env python3

"""
Convert OSM GeoPackage buildings to STL for OpenFOAM.

Height priority:
    1. OSM tag: height
    2. OSM tag: building:levels * 3.2
    3. Default = 10 m

Usage:
    python gpkg2stl.py LdW.gpkg LdW.stl
"""

import re
import sys

import geopandas as gpd
import trimesh
import pyvista as pv

DEFAULT_HEIGHT_M = 10.0
LEVEL_HEIGHT_M = 3.2

def parse_height(tags):
    """
    Parse OSM height tags from other_tags.
    """

    if tags is None:
        return None

    # height=>"18"
    m = re.search(r'"height"=>"([^"]+)"', str(tags))

    if m:
        try:
            value = m.group(1).replace("m", "").strip()
            return float(value)
        except Exception:
            pass

    return None


def parse_levels(tags):
    """
    Parse OSM building:levels tags from other_tags.
    """

    if tags is None:
        return None

    m = re.search(
        r'"building:levels"=>"([^"]+)"',
        str(tags)
    )

    if m:
        try:
            return float(m.group(1))
        except Exception:
            pass

    return None


def determine_height(row):
    """
    Determine building height.
    """

    tags = row.get("other_tags")

    h = parse_height(tags)

    if h is not None:
        return h, "height"

    levels = parse_levels(tags)

    if levels is not None:
        return levels * LEVEL_HEIGHT_M, "levels"

    return DEFAULT_HEIGHT_M, "default"


def extract_buildings(gpkg_file):

    gdf = gpd.read_file(
        gpkg_file,
        layer="multipolygons"
    )

    buildings = gdf[
        gdf["building"].notna()
    ].copy()

    if len(buildings) == 0:
        raise RuntimeError(
            "No buildings found."
        )

    print(
        f"Found {len(buildings)} building footprints"
    )

    return buildings


def project_to_metric(buildings):

    print(
        "Reprojecting to EPSG:25832 ..."
    )

    return buildings.to_crs(25832)


def build_mesh(buildings):

    meshes = []

    n_height = 0
    n_levels = 0
    n_default = 0
    n_failed = 0

    for _, row in buildings.iterrows():

        height, source = determine_height(row)

        if source == "height":
            n_height += 1
        elif source == "levels":
            n_levels += 1
        else:
            n_default += 1

        try:

            geom = row.geometry

            if geom.geom_type == "Polygon":

                mesh = trimesh.creation.extrude_polygon( geom, height, engine="earcut" )
                
                meshes.append(mesh)

            elif geom.geom_type == "MultiPolygon":

                for poly in geom.geoms:

                    mesh = trimesh.creation.extrude_polygon(
                        poly,
                        height
                    )

                    meshes.append(mesh)

        except Exception:
            n_failed += 1
            print( f"Failed geometry {_}: " f"{type(e).__name__}: {e}")

    print("\nHeight statistics:")
    print(
        f"  from 'height'          : {n_height}"
    )
    print(
        f"  from 'building:levels' : {n_levels}"
    )
    print(
        f"  default {DEFAULT_HEIGHT_M:.1f} m      : {n_default}"
    )
    print(
        f"  failed geometries      : {n_failed}"
    )

    if not meshes:
        raise RuntimeError(
            "No geometry generated."
        )

    return trimesh.util.concatenate(meshes)


def main():

    if len(sys.argv) != 3:

        print(
            "Usage:\n"
            "  python gpkg2stl.py input.gpkg output.stl"
        )

        sys.exit(1)

    input_gpkg = sys.argv[1]
    output_stl = sys.argv[2]

    buildings = extract_buildings(input_gpkg)
    buildings["geometry"] = ( buildings.geometry.buffer(0) )

    buildings = project_to_metric(buildings)
    buildings["geometry"] = ( buildings.geometry.buffer(0) )

    city = build_mesh(buildings)

    print( f"\nVertices: {len(city.vertices)}" )

    print( f"Faces   : {len(city.faces)}" )

    city.export(output_stl)

    print( f"\nWritten: {output_stl}" )

    stl_file = output_stl
    png_file = output_stl.replace(".stl", "IsoView.png")

    mesh = pv.read(stl_file)

    plotter = pv.Plotter(off_screen=True)
    plotter.add_mesh( mesh, color="lightgray", show_edges=False )
    plotter.add_axes()
    plotter.camera_position = "iso"
    plotter.show(  screenshot=png_file    )
    print(f"Saved {png_file}")

if __name__ == "__main__":
    main()