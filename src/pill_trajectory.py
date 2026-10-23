#!/usr/bin/env python3
"""Summarize passive-tracer tracks written by the pitzDaily streamlines function."""

import argparse
import csv
import math
import re
import sys
from pathlib import Path


def _read_values(lines, start, count, convert):
    values = []
    index = start
    while len(values) < count:
        if index >= len(lines):
            raise ValueError("Unexpected end of VTK file while reading data")
        values.extend(convert(token) for token in lines[index].split())
        index += 1
    if len(values) != count:
        raise ValueError("VTK data array contains more values than expected")
    return values, index


def read_legacy_vtk(path):
    lines = path.read_text(encoding="utf-8").splitlines()

    points = None
    tracks = None
    fields = {}
    point_data_count = None
    index = 0
    while index < len(lines):
        tokens = lines[index].split()
        if not tokens:
            index += 1
            continue

        if tokens[0] == "POINTS":
            count = int(tokens[1])
            values, index = _read_values(lines, index + 1, count * 3, float)
            points = [tuple(values[i : i + 3]) for i in range(0, len(values), 3)]
            continue

        if tokens[0] == "LINES":
            count = int(tokens[1])
            total = int(tokens[2])
            values, index = _read_values(lines, index + 1, total, int)
            tracks = []
            offset = 0
            for _ in range(count):
                size = values[offset]
                offset += 1
                tracks.append(values[offset : offset + size])
                offset += size
            if offset != total:
                raise ValueError("VTK polyline connectivity is malformed")
            continue

        if tokens[0] == "POINT_DATA":
            point_data_count = int(tokens[1])
        elif tokens[0] == "SCALARS":
            if point_data_count is None:
                raise ValueError("VTK scalar data appears before POINT_DATA")
            name = tokens[1]
            components = int(tokens[3]) if len(tokens) > 3 else 1
            index += 1
            if index >= len(lines) or not lines[index].startswith("LOOKUP_TABLE"):
                raise ValueError("VTK SCALARS section is missing LOOKUP_TABLE")
            values, index = _read_values(
                lines, index + 1, point_data_count * components, float
            )
            fields[name] = [
                tuple(values[i : i + components])
                for i in range(0, len(values), components)
            ]
            continue
        elif tokens[0] in ("VECTORS", "NORMALS"):
            if point_data_count is None:
                raise ValueError("VTK vector data appears before POINT_DATA")
            name = tokens[1]
            values, index = _read_values(lines, index + 1, point_data_count * 3, float)
            fields[name] = [tuple(values[i : i + 3]) for i in range(0, len(values), 3)]
            continue

        index += 1

    if points is None or tracks is None:
        raise ValueError("VTK file must contain POINTS and LINES sections")
    if point_data_count != len(points):
        raise ValueError("VTK POINT_DATA count does not match POINTS")
    if "age" not in fields or "U" not in fields:
        raise ValueError("VTK file must contain point data arrays named age and U")

    return points, tracks, fields


def outlet_x_from_block_mesh(case_dir):
    block_mesh = (case_dir / "system" / "blockMeshDict").read_text(
        encoding="utf-8"
    )
    block_mesh = re.sub(r"/\*.*?\*/|//[^\n]*", "", block_mesh, flags=re.S)
    scale_match = re.search(r"\bconvertToMeters\s+([+\-0-9.eE]+)\s*;", block_mesh)
    vertices_match = re.search(r"\bvertices\s*\((.*?)\)\s*;", block_mesh, re.S)
    if scale_match is None or vertices_match is None:
        raise ValueError("Could not read convertToMeters and vertices from blockMeshDict")

    vertices = re.findall(
        r"\(\s*([+\-0-9.eE]+)\s+([+\-0-9.eE]+)\s+([+\-0-9.eE]+)\s*\)",
        vertices_match.group(1),
    )
    if not vertices:
        raise ValueError("No vertices found in blockMeshDict")
    return max(float(vertex[0]) for vertex in vertices) * float(scale_match.group(1))


def evaluate_tracks(points, tracks, fields, outlet_x, max_time, stuck_speed):
    outcomes = []
    for track_id, point_ids in enumerate(tracks):
        if not point_ids:
            continue

        track_points = [points[point_id] for point_id in point_ids]
        ages = [fields["age"][point_id][0] for point_id in point_ids]
        velocities = [fields["U"][point_id] for point_id in point_ids]
        within_limit = [
            i for i, age in enumerate(ages) if math.isfinite(age) and age <= max_time
        ]
        if not within_limit:
            within_limit = [0]

        end_index = within_limit[-1]
        min_speed = min(
            math.sqrt(sum(component * component for component in velocities[i]))
            for i in within_limit
        )
        reached_indices = [
            i
            for i in within_limit
            if track_points[i][0] >= outlet_x - 1e-7
        ]
        reached = bool(reached_indices)
        elapsed = ages[reached_indices[0]] if reached else ages[end_index]
        outcome = (
            "reached_outlet"
            if reached
            else "stuck_low_speed"
            if min_speed <= stuck_speed
            else "not_reached"
        )
        start = track_points[0]
        end = track_points[reached_indices[0] if reached else end_index]

        outcomes.append(
            {
                "track_id": track_id,
                "outcome": outcome,
                "elapsed_time_s": elapsed,
                "min_speed_m_s": min_speed,
                "start_x_m": start[0],
                "start_y_m": start[1],
                "start_z_m": start[2],
                "end_x_m": end[0],
                "end_y_m": end[1],
                "end_z_m": end[2],
            }
        )
    return outcomes


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Analyze passive-tracer streamlines from the pitzDaily inlet. "
            "Run Allrun and postProcess -func pillTrajectories first."
        )
    )
    parser.add_argument("--case", type=Path, required=True, help="pitzDaily case directory")
    parser.add_argument(
        "--residence-time",
        type=float,
        required=True,
        help="Fluid residence time in the domain, in seconds",
    )
    parser.add_argument(
        "--stuck-speed",
        type=float,
        default=0.01,
        help="Low-speed threshold for stuck classification, in m/s (default: 0.01)",
    )
    parser.add_argument(
        "--tracks",
        type=Path,
        help="VTK tracks file (default: latest pillTrajectories output)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="CSV output path (default: case/postProcessing/pillTrajectoryResults.csv)",
    )
    args = parser.parse_args(argv)

    if not math.isfinite(args.residence_time) or args.residence_time <= 0:
        parser.error("--residence-time must be a finite positive number")
    if not math.isfinite(args.stuck_speed) or args.stuck_speed < 0:
        parser.error("--stuck-speed must be finite and non-negative")

    case_dir = args.case.resolve()
    if args.tracks:
        tracks_path = args.tracks.resolve()
    else:
        track_files = list(
            (case_dir / "postProcessing" / "pillTrajectories").glob("*/tracks.vtk")
        )
        if not track_files:
            parser.error(
                "No tracks.vtk found; run postProcess -func pillTrajectories "
                "after solving the case"
            )
        tracks_path = max(
            track_files,
            key=lambda path: float(path.parent.name),
        )

    try:
        points, tracks, fields = read_legacy_vtk(tracks_path)
        outcomes = evaluate_tracks(
            points,
            tracks,
            fields,
            outlet_x_from_block_mesh(case_dir),
            100 * args.residence_time,
            args.stuck_speed,
        )
        if not outcomes:
            raise ValueError("The VTK file contains no trajectory tracks")

        output_path = args.output or (
            case_dir / "postProcessing" / "pillTrajectoryResults.csv"
        )
        if not output_path.is_absolute():
            output_path = case_dir / output_path
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", newline="", encoding="utf-8") as output_file:
            writer = csv.DictWriter(output_file, fieldnames=outcomes[0].keys())
            writer.writeheader()
            writer.writerows(outcomes)

    except (OSError, ValueError, IndexError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1

    reached = [row for row in outcomes if row["outcome"] == "reached_outlet"]
    stuck = [row for row in outcomes if row["outcome"] == "stuck_low_speed"]
    fastest = min(reached, key=lambda row: row["elapsed_time_s"], default=None)
    critical = max(stuck, key=lambda row: row["elapsed_time_s"], default=None)

    print(f"Analyzed {len(outcomes)} inlet trajectories")
    print(f"Time limit: {100 * args.residence_time:g} s")
    print(f"Reached outlet: {len(reached)}; stuck: {len(stuck)}")
    if fastest:
        print(
            "Fastest outlet trajectory: "
            f"#{fastest['track_id']} in {fastest['elapsed_time_s']:.6g} s "
            f"from y={fastest['start_y_m']:.6g} m"
        )
    else:
        print("No trajectory reached the outlet within the time limit.")
    if critical:
        print(
            "Longest-duration stuck trajectory: "
            f"#{critical['track_id']} at y={critical['start_y_m']:.6g} m "
            f"(elapsed {critical['elapsed_time_s']:.6g} s, "
            f"minimum speed {critical['min_speed_m_s']:.6g} m/s)"
        )
    else:
        print("No trajectory met the low-speed stuck criterion.")
    print(f"Per-trajectory results: {output_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
