#include "argList.H"
#include "fvCFD.H"
#include "streamlines.H"
#include "instantList.H"
#include "meshSearch.H"
#include "StreamlineSummaryWriter.H"
#include "StreamlineTrackData.H"
#include "StreamlineVtkReader.H"

using namespace Foam;

namespace
{

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
    const fileName outputBase(analysisCase/"postProcessing");
    if (!isDir(outputBase) && !mkDir(outputBase))
    {
        FatalErrorInFunction
            << "Could not create output directory " << outputBase
            << exit(FatalError);
    }
        
    //too detailed const fileName finalOutput(outputBase/targetCase.name()/timeDirs[latestTimeIndex].name();
    const fileName finalOutput(outputBase/targetCase.name());

    if (!isDir(finalOutput) && !mkDir(finalOutput)){
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
            tracks = StreamlineVtkReader::read(temporaryVtk);
        }

        meshSearch search(mesh);
        StreamlineSummaryWriter summaryWriter
        (
            search,
            mesh,
            destinationPatches
        );
        summaryWriter.write(finalOutput/"streamlines.csv", tracks);
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
