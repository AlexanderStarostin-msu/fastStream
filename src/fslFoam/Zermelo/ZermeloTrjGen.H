#ifndef ZermeloTrjGen_H
#define ZermeloTrjGen_H

#include <boost/geometry.hpp>

#include "ZermeloTrajectory.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "wordList.H"

#include <vector>

namespace Foam
{

class ZermeloTrjGen
{
public:
    typedef boost::geometry::model::d2::point_xy<scalar> Point2;
    typedef boost::geometry::model::segment<Point2> Segment2;
    typedef boost::geometry::model::polygon<Point2> Polygon2;

protected:
    const fvMesh& mesh_;
    scalar routeHeight_;
    scalar minObstacleClearance_;
    scalar deviceSpeed_;
    scalar geometryTolerance_;
    point startPoint_;
    point endPoint_;
    std::vector<Polygon2> obstacleFootprints_;
    std::vector<Polygon2> bufferedObstacles_;

    virtual bool buildTrajectory
    (
        const point& start,
        const point& end,
        std::vector<point>& waypoints
    ) const = 0;

    bool pointIsClear(const Point2& position) const;
    bool segmentIsClear(const Point2& start, const Point2& end) const;

private:
    static std::vector<Segment2> intersectPatch
    (
        const fvMesh& mesh,
        const wordList& patchNames,
        const scalar height,
        const scalar tolerance
    );

    static std::vector<Polygon2> extractObstacleFootprints
    (
        const fvMesh& mesh,
        const wordList& patchNames,
        const scalar height,
        const scalar tolerance
    );

    void prepareObstacles();

public:
    ZermeloTrjGen
    (
        const fvMesh& mesh,
        const dictionary& inputDict,
        const dictionary& generatorDict
    );

    virtual ~ZermeloTrjGen();

    void generate
    (
        const fileName& tracksVtk,
        const fileName& linesVtk
    );
};

}

#endif
