ZermeloTrjGen
 ├─ parse configuration
 ├─ validate
 ├─ intersect patches
 ├─ determine one start and destination
 ├─ prepare obstacles
 ├─ buildTrajectory() <-- virtual
 ├─ compute length
 ├─ compute travel time
 └─ collect results

ZermeloTrajectoryWriter
 ├─ write trajectory geometry
 └─ write trajectory identifiers and metadata

ZermeloTrjNoFlow
 └─ buildTrajectory()
      ├─ build visibility graph
      ├─ add obstacle vertices
      ├─ A* shortest-path search
      ├─ reconstruct path
      └─ return baseline knots

ZermeloTrjFull
 └─ buildTrajectory()
      ├─ verify/generate baseline
      ├─ extract baseline knots
      ├─ read velocity field U
      ├─ build local flow-aware graphs
      ├─ solve minimum-time segments
      ├─ preserve knot ordering
      ├─ merge segment trajectories
      ├─ smooth curve at knots
      ├─ optimize travel time
      ├─ validate obstacle clearance
      └─ return wind-corrected trajectory