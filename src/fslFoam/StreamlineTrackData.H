#ifndef StreamlineTrackData_H
#define StreamlineTrackData_H

#include "label.H"
#include "point.H"
#include "scalar.H"

#include <vector>

namespace Foam
{

struct Track
{
    std::vector<label> pointIds;
};

struct TrackData
{
    std::vector<point> points;
    std::vector<scalar> ages;
    std::vector<Track> tracks;
};

}

#endif
