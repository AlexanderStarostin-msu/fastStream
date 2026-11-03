#include "ZermeloTrajectory.H"

namespace Foam
{

ZermeloTrajectory::ZermeloTrajectory()
:
    trajectoryId(-1),
    start(point(0, 0, 0)),
    end(point(0, 0, 0)),
    waypoints(),
    length(0),
    travelTime(0)
{}

}
