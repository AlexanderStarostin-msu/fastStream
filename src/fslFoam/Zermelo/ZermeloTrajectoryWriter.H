#ifndef ZermeloTrajectoryWriter_H
#define ZermeloTrajectoryWriter_H

#include "ZermeloTrajectory.H"
#include "fileName.H"

namespace Foam
{

class ZermeloTrajectoryWriter
{
public:
    static void write
    (
        const fileName& tracksVtk,
        const fileName& linesVtk,
        const std::vector<ZermeloTrajectory>& trajectories
    );
};

}

#endif
