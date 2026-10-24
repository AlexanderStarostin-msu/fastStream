#include "argList.H"
#include "fvCFD.H"
#include "streamlines.H"
#include "instantList.H"
#include "meshSearch.H"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

using namespace Foam;

namespace
{

struct Track
{
    std::vector<label> pointIds;
};

struct TrackData
{
    std::vector<point> points;
    std::vector<scalar> ages;
    std::vector<Track> tracks;
};


void validatePatches
(
    const wordList& patchNames,
    const polyBoundaryMesh& boundary,
    const word& regionType
)
{
    if (patchNames.empty())
    {
        FatalErrorInFunction
            << "At least one " << regionType << " patch must be specified"
            << exit(FatalError);
    }

    forAll(patchNames, i)
    {
        if (boundary.findPatchID(patchNames[i]) < 0)
        {
            FatalErrorInFunction
                << "Unknown " << regionType << " patch "
                << patchNames[i] << exit(FatalError);
        }
    }
}


fileName absolutePath(const fileName& path)
{
    fileName result(path.isAbsolute() ? path : cwd()/path);
    result.clean();
    return result;
}


void createLink(const fileName& source, const fileName& destination)
{
    if (!ln(source, destination))
    {
        FatalErrorInFunction
            << "Could not create temporary case link from " << source
            << " to " << destination << exit(FatalError);
    }
}


void removeTemporaryCase
(
    const fileName& temporaryCase,
    const word& timeName
)
{
    const fileName generatedOutput(temporaryCase/"postProcessing");
    if (isDir(generatedOutput) && !rmDir(generatedOutput))
    {
        FatalErrorInFunction
            << "Could not remove temporary tracking output "
            << generatedOutput << exit(FatalError);
    }

    fileNameList links(3);
    links[0] = temporaryCase/"system";
    links[1] = temporaryCase/"constant";
    links[2] = temporaryCase/timeName;

    forAll(links, i)
    {
        if (exists(links[i]) && !rm(links[i]))
        {
            FatalErrorInFunction
                << "Could not remove temporary case link " << links[i]
                << exit(FatalError);
        }
    }

    if (!rmDir(temporaryCase))
    {
        FatalErrorInFunction
                << "Could not remove temporary case " << temporaryCase
                << exit(FatalError);
    }
}


TrackData readTracks(const fileName& vtkFile)
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


void copyFile(const fileName& source, const fileName& destination)
{
    std::ifstream input(source.c_str(), std::ios::binary);
    std::ofstream output(destination.c_str(), std::ios::binary);
    if (!input || !output)
    {
        FatalErrorInFunction
            << "Could not copy " << source << " to " << destination
            << exit(FatalError);
    }

    output << input.rdbuf();
    if (input.bad() || !output)
    {
        FatalErrorInFunction
            << "Failed while copying " << source << " to " << destination
            << exit(FatalError);
    }
}


void writeCsv
(
    const fileName& csvFile,
    const TrackData& tracks,
    const meshSearch& search,
    const fvMesh& mesh,
    const wordList& destinationPatches
)
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
        << std::setprecision(17);

    const scalar boundaryTolerance =
        max(SMALL, 1e-7*mesh.bounds().mag());

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
        const label boundaryFace = search.findNearestBoundaryFace(end);
        if (boundaryFace >= mesh.nInternalFaces())
        {
            const face& face = mesh.faces()[boundaryFace];
            const pointHit nearestPoint =
                face.nearestPoint(end, mesh.points());
            if (nearestPoint.distance() <= boundaryTolerance)
            {
                const label patchi =
                    mesh.boundaryMesh().whichPatch(boundaryFace);
                if
                (
                    patchi >= 0
                 && findIndex
                    (
                        destinationPatches,
                        mesh.boundaryMesh()[patchi].name()
                    ) >= 0
                )
                {
                    destinationPatch = mesh.boundaryMesh()[patchi].name();
                }
            }
        }

        output
            << tracki << ','
            << start.x() << ',' << start.y() << ',' << start.z() << ','
            << end.x() << ',' << end.y() << ',' << end.z() << ','
            << flightTime << ','
            << (destinationPatch.empty() ? "false" : "true") << ','
            << destinationPatch.c_str() << '\n';
    }

    if (!output)
    {
        FatalErrorInFunction
            << "Failed while writing streamline summary " << csvFile
            << exit(FatalError);
    }
}

}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Trace streamlines through the latest velocity field in an OpenFOAM case"
    );
    argList::noParallel();

    #include "setRootCase.H"
    #include "createTime.H"

    IOdictionary inputDict
    (
        IOobject
        (
            "fsStreamLineDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    string casePath;
    inputDict.lookup("case") >> casePath;
    fileName targetCase(casePath);
    if (!targetCase.isAbsolute())
    {
        targetCase = runTime.path()/targetCase;
    }
    targetCase = absolutePath(targetCase);

    if (!isDir(targetCase/"system") || !isDir(targetCase/"constant"))
    {
        FatalErrorInFunction
            << "The configured OpenFOAM case is incomplete: " << targetCase
            << exit(FatalError);
    }

    Time sourceTime
    (
        word(Time::controlDictName),
        targetCase.path(),
        fileName(targetCase.name()),
        word("system"),
        word("constant")
    );

    const instantList timeDirs(sourceTime.times());
    label latestTimeIndex = -1;
    scalar latestTimeValue = -GREAT;
    forAll(timeDirs, timei)
    {
        if
        (
            timeDirs[timei].name() != sourceTime.constant()
         && timeDirs[timei].value() > latestTimeValue
        )
        {
            latestTimeIndex = timei;
            latestTimeValue = timeDirs[timei].value();
        }
    }
    if (latestTimeIndex < 0)
    {
        FatalErrorInFunction
            << "No numeric time directories found in " << targetCase
            << exit(FatalError);
    }

    const instant latestTime(timeDirs[latestTimeIndex]);
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
        inputDict.lookupOrDefault<wordList>
        (
            "destinationPatches",
            wordList(1, word("outlet"))
        )
    );
    const label nPoints =
        inputDict.lookupOrDefault<label>("nPoints", 1000);
    const label lifeTime =
        inputDict.lookupOrDefault<label>("lifeTime", 10000);
    const label nSubCycle =
        inputDict.lookupOrDefault<label>("nSubCycle", 5);
    const Switch outputStreamLines =
        inputDict.lookupOrDefault<Switch>("outputStreamLines", false);

    if (nPoints < 1 || lifeTime < 1 || nSubCycle < 1)
    {
        FatalIOErrorInFunction(inputDict)
            << "nPoints, lifeTime, and nSubCycle must be positive"
            << exit(FatalIOError);
    }

    const fileName analysisCase(absolutePath(runTime.path()));
    const fileName outputBase(analysisCase/"postProcessing"/"fsFoam");
    if (!isDir(outputBase) && !mkDir(outputBase))
    {
        FatalErrorInFunction
            << "Could not create output directory " << outputBase
            << exit(FatalError);
    }
    const fileName finalOutput
    (
        outputBase/targetCase.name()/timeDirs[latestTimeIndex].name()
    );
    if (!isDir(finalOutput) && !mkDir(finalOutput))
    {
        FatalErrorInFunction
            << "Could not create output directory " << finalOutput
            << exit(FatalError);
    }

    fileName temporaryCase;
    for (label suffix = 0; ; ++suffix)
    {
        OStringStream name;
        name << ".tracking-" << pid() << '-' << suffix;
        temporaryCase = outputBase/name.str();
        if (!exists(temporaryCase))
        {
            break;
        }
    }
    if (!mkDir(temporaryCase))
    {
        FatalErrorInFunction
            << "Could not create temporary tracking case " << temporaryCase
            << exit(FatalError);
    }

    createLink(targetCase/"system", temporaryCase/"system");
    createLink(targetCase/"constant", temporaryCase/"constant");
    createLink
    (
        targetCase/timeDirs[latestTimeIndex].name(),
        temporaryCase/timeDirs[latestTimeIndex].name()
    );

    const fileName temporaryVtk
    (
        temporaryCase/"postProcessing"/"fsStreamlines"
       /latestTime.name()/"tracks.vtk"
    );

    TrackData tracks;
    {
        Time caseTime
        (
            word(Time::controlDictName),
            temporaryCase.path(),
            fileName(temporaryCase.name()),
            word("system"),
            word("constant")
        );
        caseTime.setTime(latestTime, latestTimeIndex);

        Info<< "Using case " << targetCase << " at time "
            << caseTime.timeName() << nl
            << "Start patches: " << startPatches << nl
            << "Destination patches: " << destinationPatches << nl << endl;

        fvMesh mesh
        (
            IOobject
            (
                fvMesh::defaultRegion,
                caseTime.timeName(),
                caseTime,
                IOobject::MUST_READ
            )
        );

        volVectorField U
        (
            IOobject
            (
                "U",
                caseTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        validatePatches(startPatches, mesh.boundaryMesh(), "start");
        validatePatches
        (
            destinationPatches,
            mesh.boundaryMesh(),
            "destination"
        );

        dictionary seedSampleSet;
        seedSampleSet.add("type", word("boundaryRandom"));
        seedSampleSet.add("patches", startPatches);
        seedSampleSet.add("nPoints", nPoints);
        seedSampleSet.add("axis", word("xyz"));

        dictionary streamlinesDict;
        streamlinesDict.add("fields", wordList(1, word("U")));
        streamlinesDict.add("U", word("U"));
        streamlinesDict.add("direction", word("forward"));
        streamlinesDict.add("lifeTime", lifeTime);
        streamlinesDict.add("nSubCycle", nSubCycle);
        streamlinesDict.add("setFormat", word("vtk"));
        streamlinesDict.add("writeFormat", word("ascii"));
        streamlinesDict.add("seedSampleSet", seedSampleSet);

        functionObjects::streamlines streamlines
        (
            "fsStreamlines",
            caseTime,
            streamlinesDict
        );
        streamlines.write();

        if (exists(temporaryVtk))
        {
            tracks = readTracks(temporaryVtk);
        }

        meshSearch search(mesh);
        writeCsv
        (
            finalOutput/"streamlines.csv",
            tracks,
            search,
            mesh,
            destinationPatches
        );
    }

    const fileName csvFile(finalOutput/"streamlines.csv");
    if (outputStreamLines && exists(temporaryVtk))
    {
        copyFile(temporaryVtk, finalOutput/"tracks.vtk");
    }

    removeTemporaryCase
    (
        temporaryCase,
        timeDirs[latestTimeIndex].name()
    );

    Info<< "Wrote streamline summary to " << csvFile << nl;
    if (outputStreamLines && exists(finalOutput/"tracks.vtk"))
    {
        Info<< "Wrote streamline VTK to " << finalOutput/"tracks.vtk" << nl;
    }
    else if (outputStreamLines)
    {
        Info<< "No streamline VTK was produced (no track points)" << nl;
    }
    else
    {
        Info<< "VTK streamline output disabled" << nl;
    }

    return 0;
}
