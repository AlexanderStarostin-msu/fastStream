#include "StreamlineSummaryWriter.H"
#include "error.H"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <vector>

namespace Foam
{

namespace
{

struct SummaryRow
{
    std::size_t streamlineId;
    point start;
    point end;
    scalar timeOfFlight;
    word destinationPatch;
};

}


StreamlineSummaryWriter::StreamlineSummaryWriter
(
    const meshSearch& search,
    const fvMesh& mesh,
    const wordList& destinationPatches
)
:
    search_(search),
    mesh_(mesh),
    destinationPatches_(destinationPatches)
{}


void StreamlineSummaryWriter::write
(
    const fileName& csvFile,
    const TrackData& tracks
) const
{
    std::ofstream output(csvFile.c_str());
    if (!output)
    {
        FatalErrorInFunction
            << "Could not write streamline summary " << csvFile
            << exit(FatalError);
    }

    output
        << "streamline_id,start_x,start_y,start_z,"
        << "end_x,end_y,end_z,time_of_flight,reached_destination,"
        << "destination_patch\n"
        << std::fixed << std::setprecision(6);

    const scalar boundaryTolerance =
        max(SMALL, 1e-7*mesh_.bounds().mag());
    std::vector<SummaryRow> rows;

    for (std::size_t tracki = 0; tracki < tracks.tracks.size(); ++tracki)
    {
        const Track& track = tracks.tracks[tracki];
        if (track.pointIds.empty())
        {
            continue;
        }

        const label startId = track.pointIds.front();
        const label endId = track.pointIds.back();
        const point& start = tracks.points[startId];
        const point& end = tracks.points[endId];
        const scalar flightTime =
            max(scalar(0), tracks.ages[endId] - tracks.ages[startId]);

        word destinationPatch(word::null);
        const label boundaryFace = search_.findNearestBoundaryFace(end);
        if (boundaryFace >= mesh_.nInternalFaces())
        {
            const face& face = mesh_.faces()[boundaryFace];
            const pointHit nearestPoint =
                face.nearestPoint(end, mesh_.points());
            if (nearestPoint.distance() <= boundaryTolerance)
            {
                const label patchi =
                    mesh_.boundaryMesh().whichPatch(boundaryFace);
                if
                (
                    patchi >= 0
                 && findIndex
                    (
                        destinationPatches_,
                        mesh_.boundaryMesh()[patchi].name()
                    ) >= 0
                )
                {
                    destinationPatch =
                        mesh_.boundaryMesh()[patchi].name();
                }
            }
        }

        SummaryRow row;
        row.streamlineId = tracki;
        row.start = start;
        row.end = end;
        row.timeOfFlight = flightTime;
        row.destinationPatch = destinationPatch;
        rows.push_back(row);
    }

    std::stable_sort
    (
        rows.begin(),
        rows.end(),
        [](const SummaryRow& a, const SummaryRow& b)
        {
            return a.timeOfFlight < b.timeOfFlight;
        }
    );

    for (std::size_t rowi = 0; rowi < rows.size(); ++rowi)
    {
        const SummaryRow& row = rows[rowi];
        output
            << row.streamlineId << ','
            << row.start.x() << ',' << row.start.y() << ','
            << row.start.z() << ','
            << row.end.x() << ',' << row.end.y() << ','
            << row.end.z() << ','
            << row.timeOfFlight << ','
            << (row.destinationPatch.empty() ? "false" : "true") << ','
            << row.destinationPatch.c_str() << '\n';
    }

    if (!output)
    {
        FatalErrorInFunction
            << "Failed while writing streamline summary " << csvFile
            << exit(FatalError);
    }
}

}
