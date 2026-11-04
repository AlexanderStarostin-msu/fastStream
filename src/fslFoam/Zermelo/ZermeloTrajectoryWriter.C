#include "ZermeloTrajectoryWriter.H"
#include "error.H"

#include <fstream>
#include <iomanip>

namespace Foam
{

namespace
{

void writeVtk
(
    const fileName& vtkFile,
    const std::vector<ZermeloTrajectory>& trajectories,
    const bool includeMetadata
)
{
    std::ofstream vtk(vtkFile.c_str());
    if (!vtk)
    {
        FatalErrorInFunction
            << "Could not write Zermelo trajectory file " << vtkFile
            << exit(FatalError);
    }

    std::size_t totalPoints = 0;
    std::size_t totalLineSize = 0;
    for
    (
        std::size_t trajectoryi = 0;
        trajectoryi < trajectories.size();
        ++trajectoryi
    )
    {
        totalPoints += trajectories[trajectoryi].waypoints.size();
        totalLineSize += trajectories[trajectoryi].waypoints.size() + 1;
    }

    //TODO: Comment trajecory according to the naming_suffix of the computational module
    // e.g.change to "Zermelo full-flow trajectories" for ZermeloTrjFull
    vtk
        << "# vtk DataFile Version 3.0\n"
        << "Zermelo no-flow trajectories\n" 
        << "ASCII\n"
        << "DATASET POLYDATA\n"
        << "POINTS " << totalPoints << " double\n"
        << std::setprecision(12);

    for
    (
        std::size_t trajectoryi = 0;
        trajectoryi < trajectories.size();
        ++trajectoryi
    )
    {
        for
        (
            std::size_t pointi = 0;
            pointi < trajectories[trajectoryi].waypoints.size();
            ++pointi
        )
        {
            const point& position =
                trajectories[trajectoryi].waypoints[pointi];
            vtk << position.x() << ' ' << position.y() << ' '
                << position.z() << '\n';
        }
    }

    vtk << "LINES " << trajectories.size() << ' ' << totalLineSize << '\n';
    std::size_t pointOffset = 0;
    for
    (
        std::size_t trajectoryi = 0;
        trajectoryi < trajectories.size();
        ++trajectoryi
    )
    {
        vtk << trajectories[trajectoryi].waypoints.size();
        for
        (
            std::size_t pointi = 0;
            pointi < trajectories[trajectoryi].waypoints.size();
            ++pointi
        )
        {
            vtk << ' ' << pointOffset++;
        }
        vtk << '\n';
    }

    if (includeMetadata && !trajectories.empty())
    {
        vtk
            << "CELL_DATA " << trajectories.size() << '\n'
            << "SCALARS trajectory_id int 1\n"
            << "LOOKUP_TABLE default\n";
        for
        (
            std::size_t trajectoryi = 0;
            trajectoryi < trajectories.size();
            ++trajectoryi
        )
        {
            vtk << trajectories[trajectoryi].trajectoryId << '\n';
        }

        vtk
            << "SCALARS trajectory_length double 1\n"
            << "LOOKUP_TABLE default\n";
        for
        (
            std::size_t trajectoryi = 0;
            trajectoryi < trajectories.size();
            ++trajectoryi
        )
        {
            vtk << trajectories[trajectoryi].length << '\n';
        }

        vtk
            << "SCALARS travel_time double 1\n"
            << "LOOKUP_TABLE default\n";
        for
        (
            std::size_t trajectoryi = 0;
            trajectoryi < trajectories.size();
            ++trajectoryi
        )
        {
            vtk << trajectories[trajectoryi].travelTime << '\n';
        }
    }

    if (!vtk)
    {
        FatalErrorInFunction
            << "Failed while writing Zermelo trajectory file " << vtkFile
            << exit(FatalError);
    }
}

}


void ZermeloTrajectoryWriter::write
(
    const fileName& tracksVtk,
    const fileName& linesVtk,
    const std::vector<ZermeloTrajectory>& trajectories
)
{
    writeVtk(tracksVtk, trajectories, false);
    writeVtk(linesVtk, trajectories, true);
}

}
