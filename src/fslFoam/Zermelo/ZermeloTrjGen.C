#include "ZermeloTrjGen.H"
#include "ZermeloTrajectoryWriter.H"
#include "DynamicList.H"
#include "error.H"
#include "OSspecific.H"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <utility>

#include <iomanip>

namespace Foam
{

namespace bg = boost::geometry;
namespace
{

typedef std::pair<long long, long long> PointKey;
typedef std::pair<label, label> EdgeKey;

void reportProgress
(
    const char* stage,
    const std::size_t completed,
    const std::size_t total,
    label& nextPercent
)
{
    const label percent =
        total == 0
      ? 100
      : static_cast<label>(100.0*completed/total);

    if (percent >= nextPercent || completed == total)
    {
        Info<< "[Zermelo] " << stage << ": " << percent << "% ("
            << completed << '/' << total << ')' << endl;
        nextPercent = (percent/10 + 1)*10;
    }
}


point centerOfSegments
(
    const std::vector<ZermeloTrjGen::Segment2>& segments,
    const scalar height,
    const scalar tolerance
)
{
    scalar totalLength = 0;
    scalar weightedX = 0;
    scalar weightedY = 0;

    for (std::size_t i = 0; i < segments.size(); ++i)
    {
        const scalar length = bg::length(segments[i]);
        totalLength += length;
        weightedX += length*0.5
        *(
            bg::get<0, 0>(segments[i])
          + bg::get<1, 0>(segments[i])
        );
        weightedY += length*0.5
        *(
            bg::get<0, 1>(segments[i])
          + bg::get<1, 1>(segments[i])
        );
    }

    if (!(totalLength > tolerance))
    {
        FatalErrorInFunction
            << "Patch intersections have zero total length"
            << exit(FatalError);
    }

    return point(weightedX/totalLength, weightedY/totalLength, height);
}


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

template<class Ring>
void appendSegmentRingIntersections
(
    const ZermeloTrjGen::Segment2& segment,
    const Ring& ring,
    const scalar dx,
    const scalar dy,
    const scalar lengthSqr,
    std::vector<scalar>& parameters
)
{
    for (std::size_t edgei = 0; edgei + 1 < ring.size(); ++edgei)
    {
        std::vector<ZermeloTrjGen::Point2> intersections;
        bg::intersection
        (
            segment,
            ZermeloTrjGen::Segment2(ring[edgei], ring[edgei + 1]),
            intersections
        );
        for
        (
            std::size_t intersectioni = 0;
            intersectioni < intersections.size();
            ++intersectioni
        )
        {
            const scalar parameter =
                (
                    (intersections[intersectioni].x()
                        - bg::get<0, 0>(segment))*dx
                  + (intersections[intersectioni].y()
                        - bg::get<0, 1>(segment))*dy
                )/lengthSqr;
            parameters.push_back(max(scalar(0), min(scalar(1), parameter)));
        }
    }
}


bool segmentEntersPolygon
(
    const ZermeloTrjGen::Point2& start,
    const ZermeloTrjGen::Point2& end,
    const ZermeloTrjGen::Polygon2& polygon,
    const scalar tolerance
)
{
    if (bg::within(start, polygon) || bg::within(end, polygon))
    {
        return true;
    }

    const scalar dx = end.x() - start.x();
    const scalar dy = end.y() - start.y();
    const scalar lengthSqr = dx*dx + dy*dy;
    if (lengthSqr <= tolerance*tolerance)
    {
        return false;
    }

    const ZermeloTrjGen::Segment2 segment(start, end);
    std::vector<scalar> parameters(1, 0);
    parameters.push_back(1);
    appendSegmentRingIntersections
    (
        segment,
        polygon.outer(),
        dx,
        dy,
        lengthSqr,
        parameters
    );
    for
    (
        std::size_t ringi = 0;
        ringi < polygon.inners().size();
        ++ringi
    )
    {
        appendSegmentRingIntersections
        (
            segment,
            polygon.inners()[ringi],
            dx,
            dy,
            lengthSqr,
            parameters
        );
    }

    const scalar parameterTolerance = tolerance/Foam::sqrt(lengthSqr);
    std::sort(parameters.begin(), parameters.end());
    parameters.erase
    (
        std::unique
        (
            parameters.begin(),
            parameters.end(),
            [parameterTolerance](const scalar a, const scalar b)
            {
                return mag(a - b) <= parameterTolerance;
            }
        ),
        parameters.end()
    );

    for
    (
        std::size_t i = 0;
        i + 1 < parameters.size();
        ++i
    )
    {
        if (parameters[i + 1] - parameters[i] <= parameterTolerance)
        {
            continue;
        }

        const scalar parameter = 0.5*(parameters[i] + parameters[i + 1]);
        const ZermeloTrjGen::Point2 midpoint
        (
            start.x() + parameter*dx,
            start.y() + parameter*dy
        );
        if (bg::within(midpoint, polygon))
        {
            return true;
        }
    }

    return false;
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
    std::size_t totalFaces = 0;

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

        totalFaces += mesh.boundaryMesh()[patchi].size();
    }

    Info<< "[Zermelo] Intersecting patches with the route plane." << endl;
    label nextPercent = 0;
    std::size_t completedFaces = 0;
    reportProgress
    (
        "Patch-plane intersections",
        completedFaces,
        totalFaces,
        nextPercent
    );

    forAll(patchNames, patchNamei)
    {
        const label patchi =
            mesh.boundaryMesh().findPatchID(patchNames[patchNamei]);
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
            ++completedFaces;
            reportProgress
            (
                "Patch-plane intersections",
                completedFaces,
                totalFaces,
                nextPercent
            );
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

    Info<< "[Zermelo] Assembling obstacle contours from " << segments.size()
        << " intersected faces." << endl;
    std::map<PointKey, label> pointLabels;
    std::vector<Point2> graphPoints;
    std::set<EdgeKey> uniqueEdges;
    label nextSegmentPercent = 0;

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
        reportProgress
        (
            "Obstacle contour assembly",
            i + 1,
            segments.size(),
            nextSegmentPercent
        );
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
    std::size_t completedContourEdges = 0;
    label nextContourPercent = 0;
    Info<< "[Zermelo] Tracing " << edges.size()
        << " obstacle contour edges." << endl;
    reportProgress
    (
        "Obstacle contour tracing",
        completedContourEdges,
        edges.size(),
        nextContourPercent
    );

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
        ++completedContourEdges;
        reportProgress
        (
            "Obstacle contour tracing",
            completedContourEdges,
            edges.size(),
            nextContourPercent
        );

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
            ++completedContourEdges;
            reportProgress
            (
                "Obstacle contour tracing",
                completedContourEdges,
                edges.size(),
                nextContourPercent
            );
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
    minObstacleClearance_
    (
        inputDict.lookupOrDefault<scalar>("minObstacleClearance", 0)
    ),
    deviceSpeed_
    (
        generatorDict.lookupOrDefault<scalar>("deviceSpeed", 0.1)
    ),
    geometryTolerance_(max(SMALL, 1e-9*mag(mesh.bounds().span()))),
    startPoint_(point(0, 0, 0)),
    endPoint_(point(0, 0, 0)),
    obstacleFootprints_(),
    bufferedObstacles_()
{
    Info<< "[Zermelo] Validating route configuration." << endl;
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
    routeHeight_ = minHeight + 0.5*(maxHeight - minHeight);
    Info<< "[Zermelo] Fixed route plane is z=" << routeHeight_ << " m."
        << endl;

    if (!(minObstacleClearance_ >= 0) || !std::isfinite(minObstacleClearance_))
    {
        FatalIOErrorInFunction(inputDict)
            << "minObstacleClearance must be non-negative and finite"
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
    const wordList obstaclePatches
    (
        inputDict.lookupOrDefault<wordList>
        (
            "obstaclePatches",
            wordList(1, word("ground"))
        )
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

    Info<< "[Zermelo] Resolving startPoint and endPoint." << endl;
    const wordList* endpointPatches[] =
    {
        &startPatches,
        &destinationPatches
    };
    const word endpointTypes[] = {word("start"), word("destination")};
    for (label typei = 0; typei < 2; ++typei)
    {
        for
        (
            label patchi = 0;
            patchi < endpointPatches[typei]->size();
            ++patchi
        )
        {
            if
            (
                mesh_.boundaryMesh().findPatchID
                (
                    (*endpointPatches[typei])[patchi]
                ) < 0
            )
            {
                FatalIOErrorInFunction(inputDict)
                    << "Unknown " << endpointTypes[typei] << " patch "
                    << (*endpointPatches[typei])[patchi]
                    << exit(FatalIOError);
            }
        }
    }

    std::vector<Segment2> startSegments;
    if (!inputDict.found("startPoint"))
    {
        startSegments =
            intersectPatch
            (
                mesh_,
                startPatches,
                routeHeight_,
                geometryTolerance_
            );
        if (startSegments.empty())
        {
            FatalIOErrorInFunction(inputDict)
                << "startPatches do not intersect the Zermelo route plane"
                << exit(FatalIOError);
        }
    }

    std::vector<Segment2> endSegments;
    if (!inputDict.found("endPoint"))
    {
        endSegments =
            intersectPatch
            (
                mesh_,
                destinationPatches,
                routeHeight_,
                geometryTolerance_
            );
        if (endSegments.empty())
        {
            FatalIOErrorInFunction(inputDict)
                << "destinationPatches do not intersect the Zermelo route "
                << "plane"
                << exit(FatalIOError);
        }
    }

    if (inputDict.found("startPoint"))
    {
        inputDict.lookup("startPoint") >> startPoint_;
    }
    else
    {
        startPoint_ =
            centerOfSegments(startSegments, routeHeight_, geometryTolerance_);
    }
    if (inputDict.found("endPoint"))
    {
        inputDict.lookup("endPoint") >> endPoint_;
    }
    else
    {
        endPoint_ =
            centerOfSegments
            (
                endSegments,
                routeHeight_,
                geometryTolerance_
            );
    }

    if
    (
        !std::isfinite(startPoint_.x())
     || !std::isfinite(startPoint_.y())
     || !std::isfinite(startPoint_.z())
     || mag(startPoint_.z() - routeHeight_) > geometryTolerance_
    )
    {
        FatalIOErrorInFunction(inputDict)
            << "startPoint must be finite and lie in the Zermelo route plane "
            << "at z=" << routeHeight_
            << exit(FatalIOError);
    }
    if
    (
        !std::isfinite(endPoint_.x())
     || !std::isfinite(endPoint_.y())
     || !std::isfinite(endPoint_.z())
     || mag(endPoint_.z() - routeHeight_) > geometryTolerance_
    )
    {
        FatalIOErrorInFunction(inputDict)
            << "endPoint must be finite and lie in the Zermelo route "
            << "plane at z=" << routeHeight_
            << exit(FatalIOError);
    }

    Info<< "[Zermelo] Resolved endpoints: startPoint " << startPoint_
        << ", endPoint " << endPoint_ << endl
        << "[Zermelo] Extracting and buffering obstacle footprints." << endl;
    obstacleFootprints_ =
        extractObstacleFootprints
        (
            mesh_,
            obstaclePatches,
            routeHeight_,
            geometryTolerance_
        );
    prepareObstacles();
    Info<< "[Zermelo] Obstacle preparation complete: "
        << obstacleFootprints_.size() << " footprints, "
        << bufferedObstacles_.size() << " buffered obstacles." << endl;
}


ZermeloTrjGen::~ZermeloTrjGen()
{}


void ZermeloTrjGen::prepareObstacles()
{
    const bg::strategy::buffer::distance_symmetric<scalar> distanceStrategy
    (
        minObstacleClearance_
    );
    const bg::strategy::buffer::side_straight sideStrategy;
    const bg::strategy::buffer::join_round joinStrategy(16);
    const bg::strategy::buffer::end_round endStrategy(16);
    const bg::strategy::buffer::point_circle pointStrategy(16);

    Info<< "[Zermelo] Preparing clearance buffers for "
        << obstacleFootprints_.size() << " obstacle footprints." << endl;
    label nextPercent = 0;
    reportProgress
    (
        "Obstacle clearance preparation",
        0,
        obstacleFootprints_.size(),
        nextPercent
    );

    for
    (
        std::size_t polyi = 0;
        polyi < obstacleFootprints_.size();
        ++polyi
    )
    {
        if (minObstacleClearance_ == 0)
        {
            bufferedObstacles_.push_back(obstacleFootprints_[polyi]);
            reportProgress
            (
                "Obstacle clearance preparation",
                polyi + 1,
                obstacleFootprints_.size(),
                nextPercent
            );
            continue;
        }

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
        reportProgress
        (
            "Obstacle clearance preparation",
            polyi + 1,
            obstacleFootprints_.size(),
            nextPercent
        );
    }
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
    for
    (
        std::size_t polyi = 0;
        polyi < bufferedObstacles_.size();
        ++polyi
    )
    {
        if
        (
            segmentEntersPolygon
            (
                start,
                end,
                bufferedObstacles_[polyi],
                geometryTolerance_
            )
        )
        {
            return false;
        }
    }
    return true;
}


void ZermeloTrjGen::generate
(
    const fileName& tracksVtk,
    const fileName& linesVtk
)
{
    std::vector<point> waypoints;
    Info<< "[Zermelo] Route generation started (0%)." << endl;
    if (!buildTrajectory(startPoint_, endPoint_, waypoints))
    {
        if (exists(tracksVtk) && !rm(tracksVtk))
        {
            FatalErrorInFunction
                << "Could not remove stale Zermelo trajectory file "
                << tracksVtk << exit(FatalError);
        }
        if (exists(linesVtk) && !rm(linesVtk))
        {
            FatalErrorInFunction
                << "Could not remove stale Zermelo trajectory file "
                << linesVtk << exit(FatalError);
        }
        Info<< "No optimal Zermelo trajectory was found from "
            << "startPoint " << startPoint_
            << " to endPoint " << endPoint_
            << "; no trajectory output was written." << nl;
        return;
    }

    if (waypoints.empty())
    {
        FatalErrorInFunction
            << "The Zermelo trajectory strategy returned an empty trajectory"
            << exit(FatalError);
    }

    scalar trajectoryLength = 0;
    for (std::size_t pointi = 1; pointi < waypoints.size(); ++pointi)
    {
        trajectoryLength += mag(waypoints[pointi] - waypoints[pointi - 1]);
    }
    Info<< "[Zermelo] Route geometry is complete (90%); length="
        << trajectoryLength << " m. Writing trajectory files." << endl;

    ZermeloTrajectory trajectory;
    trajectory.trajectoryId = 0;
    trajectory.start = startPoint_;
    trajectory.end = endPoint_;
    trajectory.waypoints = waypoints;
    trajectory.length = trajectoryLength;
    trajectory.travelTime = trajectoryLength/deviceSpeed_;

    ZermeloTrajectoryWriter::write
    (
        tracksVtk,
        linesVtk,
        std::vector<ZermeloTrajectory>(1, trajectory)
    );
    // Calculate the direct straight-line distance between start and end
    scalar lineLength = mag(trajectory.end - trajectory.start);

    Info<< "Zermelo trajectory " << trajectory.trajectoryId << nl
        << ": startPoint " << trajectory.start
        << ", endPoint " << trajectory.end
        << ", straight-line " << lineLength
        << ", length " <<  trajectory.length
        << ", travelTime " << trajectory.travelTime << nl
        << "Wrote the optimal Zermelo trajectory to " << tracksVtk
        << " and " << linesVtk
        << ". Route generation complete (100%)." << endl;

}

}
