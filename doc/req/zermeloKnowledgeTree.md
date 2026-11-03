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

ZermeloTrjNoFlow
 └─ buildTrajectory()
        visibility graph
        shortest path

ZermeloTrajectoryWriter
 ├─ write trajectory geometry
 └─ write trajectory identifiers and metadata