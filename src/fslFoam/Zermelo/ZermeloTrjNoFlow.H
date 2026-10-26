#ifndef ZermeloTrjNoFlow_H
#define ZermeloTrjNoFlow_H

#include "ZermeloTrjGen.H"

namespace Foam
{

class ZermeloTrjNoFlow
:
    public ZermeloTrjGen
{
protected:
    virtual bool buildRoute
    (
        const point& start,
        std::vector<point>& waypoints
    ) const;

public:
    ZermeloTrjNoFlow
    (
        const fvMesh& mesh,
        const dictionary& inputDict,
        const dictionary& generatorDict
    );
};

}

#endif
