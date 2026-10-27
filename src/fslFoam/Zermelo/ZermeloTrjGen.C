#include "ZermeloTrjGen.H"
#include "ZermeloRouteWriter.H"
#include "DynamicList.H"
#include "error.H"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <utility>

namespace Foam
{

namespace bg = boost::geometry;
namespace
{

typedef std::pair<long long, long long> PointKey;
typedef std::pair<label, label> EdgeKey;


PointKey pointKey
(
    const ZermeloTrjGen::Point2& position,
    const scalar tolerance
)
{
    return PointKey
    (
        static_cast<long long>(std::llround(position.x()/tolerance)),
        static_cast<long long>(std::llround(position.y()/tolerance))
    );
}


bool samePoint
(
    const ZermeloTrjGen::Point2& a,
    const ZermeloTrjGen::Point2& b,
    const scalar tolerance
)
{
    return bg::distance(a, b) <= tolerance;
}


std::vector<ZermeloTrjGen::Point2> facePlaneIntersections
(
    const face& meshFace,
    const pointField& points,
    const scalar height,
    const scalar tolerance
)
{
    std::vector<ZermeloTrjGen::Point2> intersections;

    forAll(meshFace, edgei)
    {
        const point& a = points[meshFace[edgei]];
        const point& b = points[meshFace[(edgei + 1) % meshFace.size()]];
        const scalar da = a.z() - height;
        const scalar db = b.z() - height;

        if (mag(da) <= tolerance)
        {
            const ZermeloTrjGen::Point2 candidate(a.x(), a.y());
            bool duplicate = false;
            for
            (
                std::size_t i = 0;
                i < intersections.size();
                ++i
            )
            {
                duplicate = duplicate
                 || samePoint(intersections[i], candidate, tolerance);
            }
            if (!duplicate)
            {
                intersections.push_back(candidate);
            }
        }

        if
        (
            (da < -tolerance && db > tolerance)
         || (da > tolerance && db < -tolerance)
        )
        {
            const scalar fraction = (height - a.z())/(b.z() - a.z());
            const ZermeloTrjGen::Point2 candidate
            (
                a.x() + fraction*(b.x() - a.x()),
                a.y() + fraction*(b.y() - a.y())
            );
            bool duplicate = false;
            for
            (
                std::size_t i = 0;
                i < intersections.size();
                ++i
            )
            {
                duplicate = duplicate
                 || samePoint(intersections[i], candidate, tolerance);
            }
            if (!duplicate)
            {
                intersections.push_back(candidate);
            }
        }
    }

    if (intersections.size() < 2)
    {
        return std::vector<ZermeloTrjGen::Point2>();
    }

    scalar farthestDistance = -GREAT;
    std::pair<std::size_t, std::size_t> farthestPair(0, 1);
    for (std::size_t i = 0; i < intersections.size(); ++i)
    {
        for (std::size_t j = i + 1; j < intersections.size(); ++j)
        {
            const scalar distance =
                bg::distance(intersections[i], intersections[j]);
            if (distance > farthestDistance)
            {
                farthestDistance = distance;
                farthestPair = std::make_pair(i, j);
            }
        }
    }

    if (farthestDistance <= tolerance)
    {
        return std::vector<ZermeloTrjGen::Point2>();
    }

    return std::vector<ZermeloTrjGen::Point2>
    {
        intersections[farthestPair.first],
        intersections[farthestPair.second]
    };
}


}


std::vector<ZermeloTrjGen::Segment2> ZermeloTrjGen::intersectPatch
(
    const fvMesh& mesh,
    const wordList& patchNames,
    const scalar height,
    const scalar tolerance
)
{
    std::vector<Segment2> segments;

    forAll(patchNames, patchNamei)
    {
        const label patchi =
            mesh.boundaryMesh().findPatchID(patchNames[patchNamei]);
        if (patchi < 0)
        {
            FatalErrorInFunction
                << "Unknown Zermelo patch " << patchNames[patchNamei]
                << exit(FatalError);
        }

        const polyPatch& patch = mesh.boundaryMesh()[patchi];
        const label patchEnd = patch.start() + patch.size();
        for (label facei = patch.start(); facei < patchEnd; ++facei)
        {
            const std::vector<Point2> cut =
                facePlaneIntersections
                (
                    mesh.faces()[facei],
                    mesh.points(),
                    height,
                    tolerance
                );
            if (cut.size() == 2)
            {
                segments.push_back(Segment2(cut[0], cut[1]));
            }
        }
    }

    return segments;
}


std::vector<ZermeloTrjGen::Polygon2>
ZermeloTrjGen::extractObstacleFootprints
(
    const fvMesh& mesh,
    const wordList& patchNames,
    const scalar height,
    const scalar tolerance
)
{
    const std::vector<Segment2> segments =
        intersectPatch(mesh, patchNames, height, tolerance);
    if (segments.empty())
    {
        FatalErrorInFunction
            << "Obstacle patches do not intersect the Zermelo route plane"
            << exit(FatalError);
    }

    std::map<PointKey, label> pointLabels;
    std::vector<Point2> graphPoints;
    std::set<EdgeKey> uniqueEdges;

    for (std::size_t i = 0; i < segments.size(); ++i)
    {
        const Point2 endpoints[] =
        {
            segments[i].first,
            segments[i].second
        };
        label endpointLabels[2];
        for (label endi = 0; endi < 2; ++endi)
        {
            const PointKey key = pointKey(endpoints[endi], tolerance);
            std::map<PointKey, label>::iterator iter = pointLabels.find(key);
            if (iter == pointLabels.end())
            {
                endpointLabels[endi] = graphPoints.size();
                pointLabels.insert
                (
                    std::make_pair(key, endpointLabels[endi])
                );
                graphPoints.push_back(endpoints[endi]);
            }
            else
            {
                endpointLabels[endi] = iter->second;
            }
        }

        if (endpointLabels[0] != endpointLabels[1])
        {
            uniqueEdges.insert
            (
                std::minmax(endpointLabels[0], endpointLabels[1])
            );
        }
    }

    std::vector<EdgeKey> edges(uniqueEdges.begin(), uniqueEdges.end());
    std::vector<std::vector<label>> adjacency(graphPoints.size());
    for (std::size_t edgei = 0; edgei < edges.size(); ++edgei)
    {
        adjacency[edges[edgei].first].push_back(edgei);
        adjacency[edges[edgei].second].push_back(edgei);
    }

    std::vector<bool> visited(edges.size(), false);
    std::vector<Polygon2> polygons;

    for (std::size_t edgei = 0; edgei < edges.size(); ++edgei)
    {
        if (visited[edgei])
        {
            continue;
        }

        const label startNode = edges[edgei].first;
        label previousNode = startNode;
        label currentNode = edges[edgei].second;
        std::vector<label> loop(1, startNode);
        visited[edgei] = true;

        while (currentNode != startNode)
        {
            loop.push_back(currentNode);
            if (adjacency[currentNode].size() != 2)
            {
                FatalErrorInFunction
                    << "Obstacle patch intersection is not a set of closed "
                    << "footprints (a contour vertex has "
                    << adjacency[currentNode].size() << " segments)"
                    << exit(FatalError);
            }

            label nextEdge = -1;
            for
            (
                std::size_t adjacentEdgei = 0;
                adjacentEdgei < adjacency[currentNode].size();
                ++adjacentEdgei
            )
            {
                const label candidateEdge =
                    adjacency[currentNode][adjacentEdgei];
                if (!visited[candidateEdge])
                {
                    nextEdge = candidateEdge;
                    break;
                }
            }
            if (nextEdge < 0)
            {
                FatalErrorInFunction
                    << "Obstacle patch intersection contains an open contour"
                    << exit(FatalError);
            }

            visited[nextEdge] = true;
            const EdgeKey& edge = edges[nextEdge];
            const label nextNode =
                edge.first == currentNode ? edge.second : edge.first;
            previousNode = currentNode;
            currentNode = nextNode;
            (void)previousNode;
        }

        if (loop.size() < 3)
        {
            FatalErrorInFunction
                << "Obstacle footprint has fewer than three vertices"
                << exit(FatalError);
        }

        Polygon2 polygon;
        for (std::size_t pointi = 0; pointi < loop.size(); ++pointi)
        {
            polygon.outer().push_back(graphPoints[loop[pointi]]);
        }
        polygon.outer().push_back(polygon.outer().front());
        bg::correct(polygon);

        std::string validityReason;
        if (!bg::is_valid(polygon, validityReason))
        {
            FatalErrorInFunction
                << "Invalid obstacle footprint: " << validityReason
                << exit(FatalError);
        }
        polygons.push_back(polygon);
    }

    if (polygons.empty())
    {
        FatalErrorInFunction
            << "No closed obstacle footprints were found"
            << exit(FatalError);
    }

    return polygons;
}


ZermeloTrjGen::ZermeloTrjGen
(
    const fvMesh& mesh,
    const dictionary& inputDict,
    const dictionary& generatorDict
)
:
    mesh_(mesh),
    routeHeight_(0),
    deviceRadius_
    (
        generatorDict.lookupOrDefault<scalar>("deviceRadius", -1)
    ),
    obstacleMargin_
    (
        generatorDict.lookupOrDefault<scalar>("obstacleMargin", 0)
    ),
    deviceSpeed_
    (
        generatorDict.lookupOrDefault<scalar>("deviceSpeed", 0.1)
    ),
    geometryTolerance_(max(SMALL, 1e-9*mag(mesh.bounds().span()))),
    startSegments_(),
    destinationSegments_(),
    obstacleFootprints_(),
    bufferedObstacles_(),
    random_
    (
        generatorDict.lookupOrDefault<label>("randomSeed", 0)
    )
{
    const Switch variableHeight =
        generatorDict.lookupOrDefault<Switch>("variableHeight", false);
    if (variableHeight)
    {
        FatalIOErrorInFunction(generatorDict)
            << "variableHeight routing is not implemented"
            << exit(FatalIOError);
    }

    if (!inputDict.found("minHeight") || !inputDict.found("maxHeight"))
    {
        FatalIOErrorInFunction(inputDict)
            << "ZermeloTrjNoFlow requires finite minHeight and maxHeight"
            << exit(FatalIOError);
    }

    scalar minHeight, maxHeight;
    inputDict.lookup("minHeight") >> minHeight;
    inputDict.lookup("maxHeight") >> maxHeight;
    if
    (
        !std::isfinite(minHeight) || !std::isfinite(maxHeight)
     || minHeight <= -GREAT || minHeight >= GREAT || minHeight > maxHeight
     || maxHeight <= -GREAT || maxHeight >= GREAT  
    )
    {
        FatalIOErrorInFunction(inputDict)
            << "Zermelo requires finite minHeight <= maxHeight"
            << exit(FatalIOError);
    }
    maxHeight = max(maxHeight, minHeight+1e-6); // Ensure maxHeight >= minHeight
    routeHeight_ = 0.5*(minHeight + maxHeight);

    if (!(deviceRadius_ > 0) || !std::isfinite(deviceRadius_))
    {
        FatalIOErrorInFunction(generatorDict)
            << "deviceRadius must be positive and finite"
            << exit(FatalIOError);
    }
    if (!(obstacleMargin_ >= 0) || !std::isfinite(obstacleMargin_))
    {
        FatalIOErrorInFunction(generatorDict)
            << "obstacleMargin must be non-negative and finite"
            << exit(FatalIOError);
    }
    if (!(deviceSpeed_ > 0) || !std::isfinite(deviceSpeed_))
    {
        FatalIOErrorInFunction(generatorDict)
            << "deviceSpeed must be positive and finite"
            << exit(FatalIOError);
    }

    const wordList startPatches
    (
        inputDict.lookupOrDefault<wordList>
        (
            "startPatches",
            wordList(1, word("inlet"))
        )
    );
    const wordList destinationPatches
    (
        inputDict.lookupOrDefault<wordList>(
            "destinationPatches",
            wordList(1, word("outlet"))
        )
    );
    const wordList obstaclePatches(
        generatorDict.lookupOrDefault<wordList>("obstaclePatches",wordList(1, word("ground")))
    );

    if (startPatches.empty() || destinationPatches.empty())
    {
        FatalIOErrorInFunction(inputDict)
            << "Zermelo requires non-empty startPatches and destinationPatches"
            << exit(FatalIOError);
    }
    if (obstaclePatches.empty())
    {
        FatalIOErrorInFunction(generatorDict)
            << "obstaclePatches must not be empty"
            << exit(FatalIOError);
    }

    startSegments_ =
        intersectPatch(mesh_, startPatches, routeHeight_, geometryTolerance_);
    destinationSegments_ =
        intersectPatch
        (
            mesh_,
            destinationPatches,
            routeHeight_,
            geometryTolerance_
        );
    if (startSegments_.empty())
    {
        FatalIOErrorInFunction(inputDict)
            << "startPatches do not intersect the Zermelo route plane"
            << exit(FatalIOError);
    }
    if (destinationSegments_.empty())
    {
        FatalIOErrorInFunction(inputDict)
            << "destinationPatches do not intersect the Zermelo route plane"
            << exit(FatalIOError);
    }

    obstacleFootprints_ =
        extractObstacleFootprints
        (
            mesh_,
            obstaclePatches,
            routeHeight_,
            geometryTolerance_
        );
    prepareObstacles();
}


ZermeloTrjGen::~ZermeloTrjGen()
{}


void ZermeloTrjGen::prepareObstacles()
{
    const scalar clearance = deviceRadius_ + obstacleMargin_;
    const bg::strategy::buffer::distance_symmetric<scalar> distanceStrategy
    (
        clearance
    );
    const bg::strategy::buffer::side_straight sideStrategy;
    const bg::strategy::buffer::join_round joinStrategy(16);
    const bg::strategy::buffer::end_round endStrategy(16);
    const bg::strategy::buffer::point_circle pointStrategy(16);

    for
    (
        std::size_t polyi = 0;
        polyi < obstacleFootprints_.size();
        ++polyi
    )
    {
        bg::model::multi_polygon<Polygon2> buffered;
        bg::buffer
        (
            obstacleFootprints_[polyi],
            buffered,
            distanceStrategy,
            sideStrategy,
            joinStrategy,
            endStrategy,
            pointStrategy
        );
        for (std::size_t i = 0; i < buffered.size(); ++i)
        {
            bufferedObstacles_.push_back(buffered[i]);
        }
    }
}


std::vector<point> ZermeloTrjGen::sampleStarts(const label nPoints)
{
    std::vector<point> starts;
    starts.reserve(nPoints);
    scalar totalLength = 0;
    for
    (
        std::size_t segmenti = 0;
        segmenti < startSegments_.size();
        ++segmenti
    )
    {
        totalLength += bg::length(startSegments_[segmenti]);
    }
    if (!(totalLength > geometryTolerance_))
    {
        FatalErrorInFunction
            << "Start patch intersections have zero total length"
            << exit(FatalError);
    }

    std::uniform_real_distribution<scalar> pickSegment(0, totalLength);
    std::uniform_real_distribution<scalar> pickFraction(0, 1);
    for (label pointi = 0; pointi < nPoints; ++pointi)
    {
        scalar selection = pickSegment(random_);
        std::size_t segmenti = 0;
        for (; segmenti + 1 < startSegments_.size(); ++segmenti)
        {
            const scalar segmentLength = bg::length(startSegments_[segmenti]);
            if (selection <= segmentLength)
            {
                break;
            }
            selection -= segmentLength;
        }

        const Segment2& segment = startSegments_[segmenti];
        const scalar fraction = pickFraction(random_);
        starts.push_back
        (
            point
            (
                bg::get<0, 0>(segment)
                  + fraction*(bg::get<1, 0>(segment)-bg::get<0, 0>(segment)),
                bg::get<0, 1>(segment)
                  + fraction*(bg::get<1, 1>(segment)-bg::get<0, 1>(segment)),
                routeHeight_
            )
        );
    }
    return starts;
}


bool ZermeloTrjGen::pointIsClear(const Point2& position) const
{
    for
    (
        std::size_t polyi = 0;
        polyi < bufferedObstacles_.size();
        ++polyi
    )
    {
        if (bg::within(position, bufferedObstacles_[polyi]))
        {
            return false;
        }
    }
    return true;
}


bool ZermeloTrjGen::segmentIsClear
(
    const Point2& start,
    const Point2& end
) const
{
    const Point2 midpoint
    (
        0.5*(start.x() + end.x()),
        0.5*(start.y() + end.y())
    );

    for
    (
        std::size_t polyi = 0;
        polyi < bufferedObstacles_.size();
        ++polyi
    )
    {
        if
        (
            bg::within(start, bufferedObstacles_[polyi])
         || bg::within(end, bufferedObstacles_[polyi])
         || bg::within(midpoint, bufferedObstacles_[polyi])
        )
        {
            return false;
        }

        if (bg::intersects(Segment2(start, end), bufferedObstacles_[polyi]))
        {
            const scalar segmentLength = bg::distance(start, end);
            const scalar sampleStep = max
            (
                geometryTolerance_*10,
                deviceRadius_*0.25
            );
            const label sampleCount =
                max(label(2), label(ceil(segmentLength/sampleStep)));
            bool interiorCrossing = false;
            for (label samplei = 1; samplei < sampleCount; ++samplei)
            {
                const scalar fraction =
                    scalar(samplei)/sampleCount;
                const Point2 sample
                (
                    start.x() + fraction*(end.x() - start.x()),
                    start.y() + fraction*(end.y() - start.y())
                );
                if (bg::within(sample, bufferedObstacles_[polyi]))
                {
                    interiorCrossing = true;
                    break;
                }
            }
            if (interiorCrossing)
            {
                return false;
            }
        }
    }
    return true;
}


void ZermeloTrjGen::generate
(
    const label nPoints,
    const fileName& csvFile,
    const fileName& vtkFile
)
{
    if (nPoints < 1)
    {
        FatalErrorInFunction
            << "nPoints must be positive for Zermelo route generation"
            << exit(FatalError);
    }

    const std::vector<point> starts = sampleStarts(nPoints);
    std::vector<ZermeloRoute> routes;
    label unreachableCount = 0;

    for (std::size_t starti = 0; starti < starts.size(); ++starti)
    {
        std::vector<point> waypoints;
        if (!buildRoute(starts[starti], waypoints))
        {
            ++unreachableCount;
            Info<< "Zermelo start " << starti << " at " << starts[starti]
                << " has no feasible destination route" << nl;
            continue;
        }

        scalar routeLength = 0;
        for
        (
            std::size_t pointi = 1;
            pointi < waypoints.size();
            ++pointi
        )
        {
            routeLength += mag(waypoints[pointi] - waypoints[pointi - 1]);
        }

        ZermeloRoute route;
        route.routeId = starti;
        route.start = starts[starti];
        route.end = waypoints.back();
        route.waypoints = waypoints;
        route.length = routeLength;
        route.travelTime = routeLength/deviceSpeed_;
        routes.push_back(route);
    }

    ZermeloRouteWriter::write(csvFile, vtkFile, routes);
    Info<< "Wrote " << routes.size() << " Zermelo no-flow routes to "
        << csvFile << nl
        << "Unreachable Zermelo starts: " << unreachableCount << nl;
}

}
