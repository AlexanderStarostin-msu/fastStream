#include "ZermeloTrjNoFlow.H"
#include "error.H"

#include <algorithm>
#include <limits>
#include <queue>

namespace Foam
{

namespace bg = boost::geometry;
namespace
{

struct QueueEntry
{
    scalar distance;
    label node;
};

struct QueueEntryGreater
{
    bool operator()(const QueueEntry& a, const QueueEntry& b) const
    {
        return a.distance > b.distance;
    }
};

}


ZermeloTrjNoFlow::ZermeloTrjNoFlow
(
    const fvMesh& mesh,
    const dictionary& inputDict,
    const dictionary& generatorDict
)
:
    ZermeloTrjGen(mesh, inputDict, generatorDict)
{}


bool ZermeloTrjNoFlow::buildRoute
(
    const point& start,
    std::vector<point>& waypoints
) const
{
    const Point2 start2(start.x(), start.y());
    if (!pointIsClear(start2))
    {
        return false;
    }

    std::vector<Point2> nodes;
    nodes.push_back(start2);

    for
    (
        std::size_t polyi = 0;
        polyi < bufferedObstacles_.size();
        ++polyi
    )
    {
        const Polygon2& polygon = bufferedObstacles_[polyi];
        const Polygon2::ring_type* rings[1] = {&polygon.outer()};
        for (label ringi = 0; ringi < 1; ++ringi)
        {
            const Polygon2::ring_type& ring = *rings[ringi];
            for (std::size_t pointi = 0; pointi + 1 < ring.size(); ++pointi)
            {
                const Point2& candidate = ring[pointi];
                if (pointIsClear(candidate))
                {
                    nodes.push_back(candidate);
                }
            }
        }
    }

    const label nodeCount = nodes.size();
    std::vector<scalar> distances
    (
        nodeCount,
        std::numeric_limits<scalar>::max()
    );
    std::vector<label> previous(nodeCount, -1);
    std::vector<bool> visited(nodeCount, false);
    std::priority_queue
    <
        QueueEntry,
        std::vector<QueueEntry>,
        QueueEntryGreater
    > queue;

    distances[0] = 0;
    queue.push(QueueEntry{0, 0});
    while (!queue.empty())
    {
        const QueueEntry current = queue.top();
        queue.pop();
        if (visited[current.node])
        {
            continue;
        }
        visited[current.node] = true;

        for (label nodei = 0; nodei < nodeCount; ++nodei)
        {
            if
            (
                nodei == current.node
             || visited[nodei]
             || !segmentIsClear(nodes[current.node], nodes[nodei])
            )
            {
                continue;
            }

            const scalar candidateDistance =
                current.distance + bg::distance(nodes[current.node], nodes[nodei]);
            if (candidateDistance < distances[nodei])
            {
                distances[nodei] = candidateDistance;
                previous[nodei] = current.node;
                queue.push(QueueEntry{candidateDistance, nodei});
            }
        }
    }

    scalar bestLength = std::numeric_limits<scalar>::max();
    label bestNode = -1;
    Point2 bestDestination(0, 0);

    for (label nodei = 0; nodei < nodeCount; ++nodei)
    {
        if (!visited[nodei])
        {
            continue;
        }

        for
        (
            std::size_t segmenti = 0;
            segmenti < destinationSegments_.size();
            ++segmenti
        )
        {
            const Segment2& destination = destinationSegments_[segmenti];
            const scalar dx =
                bg::get<1, 0>(destination) - bg::get<0, 0>(destination);
            const scalar dy =
                bg::get<1, 1>(destination) - bg::get<0, 1>(destination);
            const scalar lengthSqr = dx*dx + dy*dy;
            scalar fraction =
                (
                    (nodes[nodei].x() - bg::get<0, 0>(destination))*dx
                  + (nodes[nodei].y() - bg::get<0, 1>(destination))*dy
                )/lengthSqr;
            fraction = min(scalar(1), max(scalar(0), fraction));
            const Point2 target
            (
                bg::get<0, 0>(destination) + fraction*dx,
                bg::get<0, 1>(destination) + fraction*dy
            );

            if (!segmentIsClear(nodes[nodei], target))
            {
                continue;
            }

            const scalar totalDistance =
                distances[nodei] + bg::distance(nodes[nodei], target);
            if (totalDistance < bestLength)
            {
                bestLength = totalDistance;
                bestNode = nodei;
                bestDestination = target;
            }
        }
    }

    if (bestNode < 0)
    {
        return false;
    }

    std::vector<label> reversePath;
    for (label nodei = bestNode; nodei >= 0; nodei = previous[nodei])
    {
        reversePath.push_back(nodei);
        if (nodei == 0)
        {
            break;
        }
    }
    if (reversePath.empty() || reversePath.back() != 0)
    {
        return false;
    }
    std::reverse(reversePath.begin(), reversePath.end());

    waypoints.reserve(reversePath.size() + 1);
    for (std::size_t i = 0; i < reversePath.size(); ++i)
    {
        const Point2& position = nodes[reversePath[i]];
        waypoints.push_back(point(position.x(), position.y(), routeHeight_));
    }

    const point destinationPoint
    (
        bestDestination.x(),
        bestDestination.y(),
        routeHeight_
    );
    if (mag(destinationPoint - waypoints.back()) > geometryTolerance_)
    {
        waypoints.push_back(destinationPoint);
    }
    return true;
}

}
