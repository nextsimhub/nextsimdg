/*----------------------------   timemesh.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __timemesh_H
#define __timemesh_H
/*----------------------------   timemesh.h     ---------------------------*/

#include <cassert>
#include <cstddef>

namespace Nextsim {

class TimeMesh {

public:
    size_t N;
    double tmax;
    double k; // step size

    TimeMesh()
        : N(0)
        , tmax(0)
        , k(0)
    {
    }

    TimeMesh(double K)
        : k(K)
    {
        assert(k > 0);
    }

    TimeMesh(double tm, int n)
        : N(n)
        , tmax(tm)
        , k(tm / n)
    {
        assert(k > 0);
    }

    void BasicInit(double tm, int n)
    {
        assert(tm > 0);
        assert(n > 0);
        N = n;
        tmax = tm;
        k = tmax / N;
    }
};

} // namespace Nextsim

/*----------------------------   timemesh.h     ---------------------------*/
/* end of #ifndef __timemesh_H */
#endif
/*----------------------------   timemesh.h     ---------------------------*/
