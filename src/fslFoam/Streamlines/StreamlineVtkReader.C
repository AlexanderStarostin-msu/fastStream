#include "StreamlineVtkReader.H"
#include "error.H"

#include <fstream>
#include <string>

namespace Foam
{

TrackData StreamlineVtkReader::read(const fileName& vtkFile)
{
    std::ifstream input(vtkFile.c_str());
    if (!input)
    {
        FatalErrorInFunction
            << "Could not read generated streamline file " << vtkFile
            << exit(FatalError);
    }

    TrackData data;
    std::string token;
    while (input >> token)
    {
        if (token == "POINTS")
        {
            label count;
            std::string valueType;
            input >> count >> valueType;
            data.points.resize(count);
            for (label i = 0; i < count; ++i)
            {
                input
                    >> data.points[i].x()
                    >> data.points[i].y()
                    >> data.points[i].z();
            }
        }
        else if (token == "VERTICES" || token == "LINES")
        {
            label count, totalSize;
            input >> count >> totalSize;
            for (label i = 0; i < count; ++i)
            {
                label pointCount;
                input >> pointCount;
                Track track;
                track.pointIds.resize(pointCount);
                for (label pointi = 0; pointi < pointCount; ++pointi)
                {
                    input >> track.pointIds[pointi];
                }
                data.tracks.push_back(track);
            }
        }
        else if (token == "FIELD")
        {
            std::string fieldName;
            label fieldCount;
            input >> fieldName >> fieldCount;
            for (label fieldi = 0; fieldi < fieldCount; ++fieldi)
            {
                std::string name, valueType;
                label components, tuples;
                input >> name >> components >> tuples >> valueType;

                if (name == "age" && components == 1)
                {
                    data.ages.resize(tuples);
                    for (label i = 0; i < tuples; ++i)
                    {
                        input >> data.ages[i];
                    }
                }
                else
                {
                    std::string ignored;
                    for (label i = 0; i < components*tuples; ++i)
                    {
                        input >> ignored;
                    }
                }
            }
        }
    }

    if (!input.eof() || data.ages.size() != data.points.size())
    {
        FatalErrorInFunction
            << "Invalid or incomplete ASCII VTK streamline data in "
            << vtkFile << exit(FatalError);
    }

    for (std::size_t tracki = 0; tracki < data.tracks.size(); ++tracki)
    {
        for
        (
            std::size_t pointi = 0;
            pointi < data.tracks[tracki].pointIds.size();
            ++pointi
        )
        {
            const label id = data.tracks[tracki].pointIds[pointi];
            if
            (
                id < 0
             || static_cast<std::size_t>(id) >= data.points.size()
            )
            {
                FatalErrorInFunction
                    << "Invalid point index in " << vtkFile
                    << exit(FatalError);
            }
        }
    }
    return data;
}

}
