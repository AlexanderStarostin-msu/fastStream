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
    scalar estimatedTotalCost; // g + h (Estimated total cost)
    scalar costFromStart;      // g (Actual distance traveled from start)
    label node;
};

struct QueueEntryGreater
{
    bool operator()(const QueueEntry& a, const QueueEntry& b) const
    {
        return a.estimatedTotalCost > b.estimatedTotalCost;
    }
};

// Generalized progress reporter with a customizable step increment (e.g., 1 for 1%, 10 for 10%)
void reportProgress
(
    const char* stage,
    const std::size_t completed,
    const std::size_t total,
    label& nextPercent,
    const int stepIncrement = 10
)
{
    const label percent =
        total == 0
      ? 100
      : static_cast<label>(100.0 * completed / total);

    if (percent >= nextPercent || completed == total)
    {
        Info<< "[Zermelo] " << stage << ": " << percent << "% ("
            << completed << '/' << total << ')' << endl;
        nextPercent = ((percent / stepIncrement) + 1) * stepIncrement;
    }
}


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
    Info<< "[Zermelo] Checking route endpoints against obstacles." << endl;
    const Point2 start2(start.x(), start.y());
    const Point2 end2(end.x(), end.y());
    if (!pointIsClear(start2) || !pointIsClear(end2))
    {
        Info<< "[Zermelo] Route construction stopped: an endpoint lies "
            << "inside an expanded obstacle." << endl;
        return false;
    }

    std::vector<Point2> nodes;
    nodes.push_back(start2);
    nodes.push_back(end2);

    std::size_t totalObstacleVertices = 0;
    for (std::size_t polyi = 0; polyi < bufferedObstacles_.size(); ++polyi)
    {
        totalObstacleVertices += bufferedObstacles_[polyi].outer().size() - 1;
    }
    Info<< "[Zermelo] Adding " << totalObstacleVertices
        << " obstacle vertices to the visibility graph." << endl;
    label nextVertexPercent = 0;
    reportProgress
    (
        "Visibility-graph vertex preparation",
        0,
        totalObstacleVertices,
        nextVertexPercent
    );
    std::size_t completedObstacleVertices = 0;
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
                ++completedObstacleVertices;
                reportProgress
                (
                    "Visibility-graph vertex preparation",
                    completedObstacleVertices,
                    totalObstacleVertices,
                    nextVertexPercent
                );
            }
        }
    }

    const label nodeCount = nodes.size();
    Info<< "[Zermelo] Visibility graph contains " << nodeCount
        << " candidate vertices" << nl
        << " Starting A* shortest-path search. Progress could be underestimated." << endl;
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
    const scalar initialH = bg::distance(nodes[0], end2);
    queue.push(QueueEntry{initialH, 0, 0});

    const std::size_t totalCandidateChecks =
        static_cast<std::size_t>(nodeCount)*nodeCount;
    std::size_t completedCandidateChecks = 0;
    label nextSearchPercent = 0;
    reportProgress
    (
        "Visibility search candidate checks",
        completedCandidateChecks,
        totalCandidateChecks,
        nextSearchPercent,
        1
    );
    std::size_t settledNodes = 0;
    
    while (!queue.empty())
    {
        const QueueEntry current = queue.top();
        queue.pop();
        
        if (visited[current.node])
        {
            continue;
        }
        visited[current.node] = true;
        ++settledNodes;
        if (current.node == 1)
        {
            const label searchPercent =
                totalCandidateChecks == 0
              ? 100
              : static_cast<label>
                (
                    100.0*completedCandidateChecks/totalCandidateChecks
                );
            Info<< "[Zermelo] Destination settled after "
                << settledNodes << '/' << nodeCount
                << " vertices and " << searchPercent
                << "% of the maximum candidate checks; A* shortest-path "
                << "optimality is proven." << endl;
            break;
        }

        for (label nodei = 0; nodei < nodeCount; ++nodei)
        {
            ++completedCandidateChecks;
            reportProgress
            (
                "Visibility search candidate checks",
                completedCandidateChecks,
                totalCandidateChecks,
                nextSearchPercent
            );
            
            if
            (
                nodei == current.node
             || visited[nodei]
             || !segmentIsClear(nodes[current.node], nodes[nodei])
            )
            {
                continue;
            }

            const scalar candidateCostFromStart =
                current.costFromStart + bg::distance(nodes[current.node], nodes[nodei]);
                
            if (candidateCostFromStart < distances[nodei])
            {
                distances[nodei] = candidateCostFromStart;
                previous[nodei] = current.node;
                
                // A* heuristic: add straight-line distance from candidate node to destination (end2)
                const scalar heuristic = bg::distance(nodes[nodei], end2);
                const scalar estimatedTotal = candidateCostFromStart + heuristic;
                
                queue.push(QueueEntry{estimatedTotal, candidateCostFromStart, nodei});
            }
        }
    }

    if (!visited[1])
    {
        const label searchPercent =
            totalCandidateChecks == 0
          ? 100
          : static_cast<label>
            (
                100.0*completedCandidateChecks/totalCandidateChecks
            );
        Info<< "[Zermelo] Visibility search exhausted "
            << settledNodes << " reachable vertices after " << searchPercent
            << "% of the maximum candidate checks; no route exists." << endl;
        return false;
    }

    Info<< "[Zermelo] Reconstructing the optimal route from predecessor "
        << "vertices." << endl;
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
    Info<< "[Zermelo] Route reconstruction complete (100%): "
        << waypoints.size() << " waypoints." << endl;
    return true;
}

}