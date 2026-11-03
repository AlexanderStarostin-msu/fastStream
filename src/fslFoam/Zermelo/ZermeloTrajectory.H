#ifndef ZermeloTrajectory_H
#define ZermeloTrajectory_H

#include "point.H"
#include "scalar.H"

#include <vector>

namespace Foam
{

struct ZermeloTrajectory
{
    label trajectoryId;
    point start;
    point end;
    std::vector<point> waypoints;
    scalar length;
    scalar travelTime;

    ZermeloTrajectory();
};

}

#endif
