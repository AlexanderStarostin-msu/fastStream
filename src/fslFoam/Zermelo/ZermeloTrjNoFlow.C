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


bool ZermeloTrjNoFlow::buildTrajectory
(
    const point& start,
    const point& end,
    std::vector<point>& waypoints
) const
{
    const Point2 start2(start.x(), start.y());
    const Point2 end2(end.x(), end.y());
    if (!pointIsClear(start2) || !pointIsClear(end2))
    {
        return false;
    }

    std::vector<Point2> nodes;
    nodes.push_back(start2);
    nodes.push_back(end2);

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

    if (!visited[1])
    {
        return false;
    }

    std::vector<label> reversePath;
    for (label nodei = 1; nodei >= 0; nodei = previous[nodei])
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

    waypoints.reserve(reversePath.size());
    for (std::size_t i = 0; i < reversePath.size(); ++i)
    {
        const Point2& position = nodes[reversePath[i]];
        waypoints.push_back(point(position.x(), position.y(), routeHeight_));
    }
    return true;
}

}
