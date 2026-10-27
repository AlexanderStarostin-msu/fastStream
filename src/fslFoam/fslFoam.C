#include "Zermelo/ZermeloTrjNoFlow.H"

#include "argList.H"
#include "fvCFD.H"
#include "streamlines.H"
#include "instantList.H"
#include "meshSearch.H"
#include "Streamlines/MotionCorridorConstraint.H"
#include "Streamlines/StreamlineSummaryWriter.H"
#include "Streamlines/StreamlineTrackData.H"
#include "Streamlines/StreamlineVtkReader.H"

#include "Random.H"

#include <algorithm>
#include <fstream>
#include <vector>

using namespace Foam;

namespace
{

struct SeedTriangle
{
    point a;
    point b;
    point c;
    scalar cumulativeArea;
};


std::vector<point> clipPolygon
(
    const std::vector<point>& polygon,
    const scalar height,
    const bool keepAbove
)
{
    std::vector<point> clipped;
    if (polygon.empty())
    {
        return clipped;
    }

    point previous = polygon.back();
    bool previousInside =
        keepAbove ? previous.z() >= height : previous.z() <= height;
    for (std::size_t pointi = 0; pointi < polygon.size(); ++pointi)
    {
        const point& current = polygon[pointi];
        const bool currentInside =
            keepAbove ? current.z() >= height : current.z() <= height;
        if (currentInside != previousInside)
        {
            const scalar fraction =
                (height - previous.z())/(current.z() - previous.z());
            clipped.push_back(previous + fraction*(current - previous));
        }
        if (currentInside)
        {
            clipped.push_back(current);
        }
        previous = current;
        previousInside = currentInside;
    }

    return clipped;
}


pointField sampleHeightConstrainedSeeds
(
    const fvMesh& mesh,
    const wordList& patchNames,
    const label nPoints,
    const scalar minHeight,
    const scalar maxHeight,
    const bool hasMinHeight,
    const bool hasMaxHeight
)
{
    std::vector<SeedTriangle> triangles;
    scalar totalArea = 0;
    // Keep sampled points inside the owner-cell tetrahedra for tracking.
    const scalar inwardFraction = 1e-4;

    forAll(patchNames, patchNamei)
    {
        const label patchi =
            mesh.boundaryMesh().findPatchID(patchNames[patchNamei]);
        const polyPatch& patch = mesh.boundaryMesh()[patchi];
        forAll(patch, patchFacei)
        {
            const face& patchFace = patch[patchFacei];
            const label facei = patchFacei + patch.start();
            const point& cellCentre =
                mesh.C()[mesh.faceOwner()[facei]];
            for (label tetPti = 1; tetPti < patchFace.size() - 1; ++tetPti)
            {
                std::vector<point> polygon(3);
                polygon[0] =
                    (1 - inwardFraction)*mesh.points()[patchFace[0]]
                  + inwardFraction*cellCentre;
                polygon[1] =
                    (1 - inwardFraction)*mesh.points()[patchFace[tetPti]]
                  + inwardFraction*cellCentre;
                polygon[2] =
                    (1 - inwardFraction)*mesh.points()[patchFace[tetPti + 1]]
                  + inwardFraction*cellCentre;

                if (hasMinHeight)
                {
                    polygon = clipPolygon(polygon, minHeight, true);
                }
                if (hasMaxHeight)
                {
                    polygon = clipPolygon(polygon, maxHeight, false);
                }
                if (polygon.size() < 3)
                {
                    continue;
                }

                for
                (
                    std::size_t pointi = 1;
                    pointi + 1 < polygon.size();
                    ++pointi
                )
                {
                    const point& a = polygon[0];
                    const point& b = polygon[pointi];
                    const point& c = polygon[pointi + 1];
                    const scalar area = 0.5*mag((b - a) ^ (c - a));
                    if (area > VSMALL)
                    {
                        totalArea += area;
                        triangles.push_back
                        (
                            SeedTriangle{a, b, c, totalArea}
                        );
                    }
                }
            }
        }
    }

    if (triangles.empty() || !(totalArea > VSMALL))
    {
        FatalErrorInFunction
            << "The configured start patches have no area within the height "
            << "range" << exit(FatalError);
    }

    Random random(261782);
    pointField seeds(nPoints);
    forAll(seeds, pointi)
    {
        const scalar selection = random.scalar01()*totalArea;
        const std::vector<SeedTriangle>::const_iterator triangle =
            std::lower_bound
            (
                triangles.begin(),
                triangles.end(),
                selection,
                [](const SeedTriangle& candidate, const scalar value)
                {
                    return candidate.cumulativeArea < value;
                }
            );
        const SeedTriangle& selected =
            triangle == triangles.end() ? triangles.back() : *triangle;
        const scalar root = Foam::sqrt(random.scalar01());
        const scalar fraction = random.scalar01();
        seeds[pointi] =
            (1 - root)*selected.a
          + root*(1 - fraction)*selected.b
          + root*fraction*selected.c;
    }

    return seeds;
}


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
    const bool hasMinHeight = inputDict.found("minHeight");
    const bool hasMaxHeight = inputDict.found("maxHeight");
    const bool heightLimited = hasMinHeight || hasMaxHeight;

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
    const fileName finalOutput(outputBase);

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

        MotionCorridorConstraint corridorConstraint(mesh, inputDict);

        if (inputDict.found("ZermeloTrjNoFlow"))
        {
            const dictionary& zermeloDict =
                inputDict.subDict("ZermeloTrjNoFlow");
            ZermeloTrjNoFlow zermelo(mesh, inputDict, zermeloDict);
            zermelo.generate
            (
                nPoints,
                finalOutput/"zermelo-no-flow.csv",
                finalOutput/"zermelo-no-flow.vtk"
            );
        }

        dictionary seedSampleSet;
        if (heightLimited)
        {
            const pointField seeds =
                sampleHeightConstrainedSeeds
                (
                    mesh,
                    startPatches,
                    nPoints,
                    inputDict.lookupOrDefault<scalar>("minHeight", -GREAT),
                    inputDict.lookupOrDefault<scalar>("maxHeight", GREAT),
                    hasMinHeight,
                    hasMaxHeight
                );
            seedSampleSet.add("type", word("points"));
            seedSampleSet.add("points", seeds);
            seedSampleSet.add("ordered", Switch(false));
            seedSampleSet.add("axis", word("xyz"));
        }
        else
        {
            seedSampleSet.add("type", word("boundaryRandom"));
            seedSampleSet.add("patches", startPatches);
            seedSampleSet.add("nPoints", nPoints);
            seedSampleSet.add("axis", word("xyz"));
        }

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
            destinationPatches,
            corridorConstraint
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
