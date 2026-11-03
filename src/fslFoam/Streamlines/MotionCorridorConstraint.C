#include "MotionCorridorConstraint.H"
#include "DynamicList.H"
#include "error.H"
#include "treeBoundBox.H"

#include <cmath>

namespace Foam
{

namespace
{

labelList obstacleFaceIds
(
    const fvMesh& mesh,
    const wordList& obstaclePatches
)
{
    DynamicList<label> faceIds;

    forAll(obstaclePatches, patchi)
    {
        const label boundaryPatchi =
            mesh.boundaryMesh().findPatchID(obstaclePatches[patchi]);
        if (boundaryPatchi < 0)
        {
            FatalErrorInFunction
                << "Unknown obstacle patch " << obstaclePatches[patchi]
                << exit(FatalError);
        }

        const polyPatch& patch = mesh.boundaryMesh()[boundaryPatchi];
        const label patchEnd = patch.start() + patch.size();
        for (label facei = patch.start(); facei < patchEnd; ++facei)
        {
            faceIds.append(facei);
        }
    }

    return faceIds;
}

}


MotionCorridorConstraint::MotionCorridorConstraint
(
    const fvMesh& mesh,
    const dictionary& inputDict
)
:
    minHeight_(inputDict.lookupOrDefault<scalar>("minHeight", -GREAT)),
    maxHeight_(inputDict.lookupOrDefault<scalar>("maxHeight", GREAT)),
    minObstacleClearance_
    (
        inputDict.lookupOrDefault<scalar>("minObstacleClearance", 0)
    ),
    obstacleFaceIds_
    (
        obstacleFaceIds
        (
            mesh,
            inputDict.lookupOrDefault<wordList>
            (
                "obstaclePatches",
                wordList()
            )
        )
    ),
    obstaclePatch_
    (
        IndirectList<face>(mesh.faces(), obstacleFaceIds_),
        mesh.points()
    ),
    obstacleTreeData_
    (
        false,
        obstaclePatch_,
        indexedOctree<obstacleTreeData>::perturbTol()
    ),
    obstacleTreePtr_(nullptr)
{
    const bool hasMinHeight = inputDict.found("minHeight");
    const bool hasMaxHeight = inputDict.found("maxHeight");

    if
    (
        hasMinHeight
     &&
        (
            !std::isfinite(minHeight_)
         || minHeight_ <= -GREAT
         || minHeight_ >= GREAT
        )
    )
    {
        FatalIOErrorInFunction(inputDict)
            << "minHeight must be finite"
            << exit(FatalIOError);
    }

    if(hasMaxHeight && ( !std::isfinite(maxHeight_) || maxHeight_ <= -GREAT || maxHeight_ >= GREAT) )
    {
        FatalIOErrorInFunction(inputDict)
            << "maxHeight must be finite"
            << exit(FatalIOError);
    }

    if (hasMinHeight && hasMaxHeight && !(minHeight_ <= maxHeight_))
    {
        FatalIOErrorInFunction(inputDict)
            << "minHeight must be less or equal to maxHeight"
            << exit(FatalIOError);
    }
    maxHeight_ = max(maxHeight_, minHeight_+1e-6); // Ensure maxHeight >= minHeight
    if (!(minObstacleClearance_ >= 0))
    {
        FatalIOErrorInFunction(inputDict)
            << "minObstacleClearance must be non-negative"
            << exit(FatalIOError);
    }

    if (minObstacleClearance_ > 0 && obstacleFaceIds_.empty())
    {
        FatalIOErrorInFunction(inputDict)
            << "minObstacleClearance requires one or more obstaclePatches"
            << exit(FatalIOError);
    }

    if (minObstacleClearance_ > 0)
    {
        const treeBoundBox overallBb(mesh.points());
        const scalar boxExtension =
            max(scalar(1e-6), scalar(1e-6)*mag(overallBb.span()));

        obstacleTreePtr_.reset
        (
            new indexedOctree<obstacleTreeData>
            (
                obstacleTreeData_,
                overallBb.extend(boxExtension),
                10,
                10,
                3.0
            )
        );
    }
}


bool MotionCorridorConstraint::isAdmissible(const point& position) const
{
    if (position.z() < minHeight_ || position.z() > maxHeight_)
    {
        return false;
    }

    if (obstacleTreePtr_.valid())
    {
        const pointIndexHit nearestObstacle = obstacleTreePtr_().findNearest
        (
            position,
            sqr(minObstacleClearance_)
        );
        if (nearestObstacle.hit())
        {
            return false;
        }
    }

    return true;
}
}
