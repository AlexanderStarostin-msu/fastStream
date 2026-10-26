#include "ZermeloRouteWriter.H"
#include "error.H"

#include <fstream>
#include <iomanip>

namespace Foam
{

void ZermeloRouteWriter::write
(
    const fileName& csvFile,
    const fileName& vtkFile,
    const std::vector<ZermeloRoute>& routes
)
{
    std::ofstream csv(csvFile.c_str());
    if (!csv)
    {
        FatalErrorInFunction
            << "Could not write Zermelo route summary " << csvFile
            << exit(FatalError);
    }

    csv
        << "route_id,start_x,start_y,start_z,end_x,end_y,end_z,"
        << "route_length,travel_time,waypoints\n"
        << std::fixed << std::setprecision(6);

    for (std::size_t routei = 0; routei < routes.size(); ++routei)
    {
        const ZermeloRoute& route = routes[routei];
        csv
            << route.routeId << ','
            << route.start.x() << ',' << route.start.y() << ','
            << route.start.z() << ','
            << route.end.x() << ',' << route.end.y() << ','
            << route.end.z() << ','
            << route.length << ',' << route.travelTime << ',';

        for
        (
            std::size_t pointi = 0;
            pointi < route.waypoints.size();
            ++pointi
        )
        {
            if (pointi)
            {
                csv << ';';
            }
            csv
                << route.waypoints[pointi].x() << ' '
                << route.waypoints[pointi].y() << ' '
                << route.waypoints[pointi].z();
        }
        csv << '\n';
    }

    if (!csv)
    {
        FatalErrorInFunction
            << "Failed while writing Zermelo route summary " << csvFile
            << exit(FatalError);
    }

    std::ofstream vtk(vtkFile.c_str());
    if (!vtk)
    {
        FatalErrorInFunction
            << "Could not write Zermelo route geometry " << vtkFile
            << exit(FatalError);
    }

    std::size_t totalPoints = 0;
    std::size_t totalLineSize = 0;
    for (std::size_t routei = 0; routei < routes.size(); ++routei)
    {
        totalPoints += routes[routei].waypoints.size();
        totalLineSize += routes[routei].waypoints.size() + 1;
    }

    vtk
        << "# vtk DataFile Version 3.0\n"
        << "Zermelo no-flow routes\n"
        << "ASCII\n"
        << "DATASET POLYDATA\n"
        << "POINTS " << totalPoints << " double\n"
        << std::setprecision(12);

    for (std::size_t routei = 0; routei < routes.size(); ++routei)
    {
        for
        (
            std::size_t pointi = 0;
            pointi < routes[routei].waypoints.size();
            ++pointi
        )
        {
            const point& position = routes[routei].waypoints[pointi];
            vtk << position.x() << ' ' << position.y() << ' '
                << position.z() << '\n';
        }
    }

    vtk << "LINES " << routes.size() << ' ' << totalLineSize << '\n';
    std::size_t pointOffset = 0;
    for (std::size_t routei = 0; routei < routes.size(); ++routei)
    {
        vtk << routes[routei].waypoints.size();
        for
        (
            std::size_t pointi = 0;
            pointi < routes[routei].waypoints.size();
            ++pointi
        )
        {
            vtk << ' ' << pointOffset++;
        }
        vtk << '\n';
    }

    if (!vtk)
    {
        FatalErrorInFunction
            << "Failed while writing Zermelo route geometry " << vtkFile
            << exit(FatalError);
    }
}

}
