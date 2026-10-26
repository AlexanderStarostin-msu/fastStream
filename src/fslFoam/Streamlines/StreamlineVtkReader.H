#ifndef StreamlineVtkReader_H
#define StreamlineVtkReader_H

#include "StreamlineTrackData.H"
#include "fileName.H"

namespace Foam
{

class StreamlineVtkReader
{
public:
    static TrackData read(const fileName& vtkFile);
};

}

#endif
