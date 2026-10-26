#ifndef ZermeloRoute_H
#define ZermeloRoute_H

#include "point.H"
#include "scalar.H"

#include <vector>

namespace Foam
{

struct ZermeloRoute
{
    label routeId;
    point start;
    point end;
    std::vector<point> waypoints;
    scalar length;
    scalar travelTime;
};

}

#endif
