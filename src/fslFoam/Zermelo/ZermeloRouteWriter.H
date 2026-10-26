#ifndef ZermeloRouteWriter_H
#define ZermeloRouteWriter_H

#include "ZermeloRoute.H"
#include "fileName.H"

namespace Foam
{

class ZermeloRouteWriter
{
public:
    static void write
    (
        const fileName& csvFile,
        const fileName& vtkFile,
        const std::vector<ZermeloRoute>& routes
    );
};

}

#endif
