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
    size_t N; //!< Number of time steps
    size_t Nsub; //!< Number of sybcycling steps in each time step
    double tmax;
    double dt; // step size
    double dt_momentum; // subcycling step size

    TimeMesh()
        : N(0)
        , tmax(0)
        , dt(0)
        , dt_momentum(0)
    {
    }

    TimeMesh(double DT, double DTmomentum)
        : dt(DT)
        , dt_momentum(DTmomentum)
    {
        assert(dt > 0);
        assert(dt_momentum > 0);
        assert(dt_momentum <= dt);

        Nsub = static_cast<size_t>(1.e-8 + dt / dt_momentum);
        assert(Nsub >= 1);
        assert(fabs(dt - Nsub * dt_momentum) < 1.e-12);
    }

    TimeMesh(double tm, size_t n, size_t nsub)
        : N(n)
        , Nsub(nsub)
        , tmax(tm)
        , dt(tm / n)
        , dt_momentum(tm / (n * nsub))
    {
        assert(dt > 0);
        assert(dt_momentum > 0);
        assert(dt_momentum < dt);
    }

    //Question
    void BasicInit(size_t n, double DT, double DTmomentum)
    {
        dt = DT;
        dt_momentum = DTmomentum;
        assert(dt > 0);
        assert(dt_momentum > 0);
        assert(dt_momentum <= dt);
        N = n;

        Nsub = static_cast<size_t>(1.e-8 + dt / dt_momentum);
        assert(Nsub >= 1);
        //assert(fabs(dt - Nsub * dt_momentum) < 1.e-12);
    }

    void BasicInit(double tm, size_t n, size_t nsub)
    {
        assert(tm > 0);
        assert(n > 0);
        assert(nsub > 0);
        N = n;
        Nsub = nsub;
        tmax = tm;
        dt = tmax / N;
        dt_momentum = tmax / (N * Nsub);
    }
};

} // namespace Nextsim

/*----------------------------   timemesh.h     ---------------------------*/
/* end of #ifndef __timemesh_H */
#endif
/*----------------------------   timemesh.h     ---------------------------*/
