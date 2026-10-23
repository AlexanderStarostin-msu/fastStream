import tempfile
import unittest
from pathlib import Path

from pill_trajectory import evaluate_tracks, outlet_x_from_block_mesh, read_legacy_vtk


class PillTrajectoryTests(unittest.TestCase):
    def test_reads_streamline_vtk_and_classifies_tracks(self):
        vtk = """# vtk DataFile Version 2.0
tracks
ASCII
DATASET POLYDATA
POINTS 4 float
0 0 0
1 0 0
0 1 0
0.5 1 0
LINES 2 6
2 0 1
2 2 3
POINT_DATA 4
SCALARS age float
LOOKUP_TABLE default
0 1 0 2
VECTORS U float
10 0 0
10 0 0
0.001 0 0
0.001 0 0
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "tracks.vtk"
            path.write_text(vtk, encoding="utf-8")
            points, tracks, fields = read_legacy_vtk(path)

        outcomes = evaluate_tracks(points, tracks, fields, 0.9, 10, 0.01)
        self.assertEqual([row["outcome"] for row in outcomes], [
            "reached_outlet",
            "stuck_low_speed",
        ])
        self.assertEqual(outcomes[0]["elapsed_time_s"], 1)
        self.assertEqual(outcomes[1]["min_speed_m_s"], 0.001)

    def test_reads_scaled_outlet_location(self):
        block_mesh = """convertToMeters 0.001;
vertices
(
    (-20.6 0 -0.5)
    (290 16.6 0.5)
);
"""
        with tempfile.TemporaryDirectory() as directory:
            case = Path(directory)
            (case / "system").mkdir()
            (case / "system" / "blockMeshDict").write_text(
                block_mesh, encoding="utf-8"
            )
            self.assertAlmostEqual(outlet_x_from_block_mesh(case), 0.29)


if __name__ == "__main__":
    unittest.main()
