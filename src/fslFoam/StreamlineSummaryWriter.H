#ifndef StreamlineSummaryWriter_H
#define StreamlineSummaryWriter_H

#include "StreamlineTrackData.H"
#include "fileName.H"
#include "fvMesh.H"
#include "meshSearch.H"
#include "wordList.H"

namespace Foam
{

class StreamlineSummaryWriter
{
    const meshSearch& search_;
    const fvMesh& mesh_;
    const wordList& destinationPatches_;

public:
    StreamlineSummaryWriter
    (
        const meshSearch& search,
        const fvMesh& mesh,
        const wordList& destinationPatches
    );

    void write(const fileName& csvFile, const TrackData& tracks) const;
};

}

#endif
