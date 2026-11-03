#ifndef MotionCorridorConstraint_H
#define MotionCorridorConstraint_H

#include "dictionary.H"
#include "fvMesh.H"
#include "indirectPrimitivePatch.H"
#include "indexedOctree.H"
#include "treeDataPrimitivePatch.H"

namespace Foam
{

class MotionCorridorConstraint
{
    typedef treeDataPrimitivePatch<indirectPrimitivePatch> obstacleTreeData;

    scalar minHeight_;
    scalar maxHeight_;
    scalar minObstacleClearance_;
    labelList obstacleFaceIds_;
    indirectPrimitivePatch obstaclePatch_;
    obstacleTreeData obstacleTreeData_;
    autoPtr<indexedOctree<obstacleTreeData>> obstacleTreePtr_;

public:
    MotionCorridorConstraint
    (
        const fvMesh& mesh,
        const dictionary& inputDict
    );

    bool isAdmissible(const point& position) const;
};

}

#endif
