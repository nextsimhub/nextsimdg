/*!
 * @file UniformOcean.cpp
 *
 * @date 23 Aug 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/UniformOcean.hpp"

#include "include/IFreezingPoint.hpp"
#include "include/NextsimModule.hpp"
#include "include/constants.hpp"

namespace Nextsim {

void UniformOcean::setData(const ModelState::DataMap& ms)
{
    IOceanBoundary::setData(ms);
    sst = sst0;
    sss = sss0;
    mld = mld0;
    u = u0;
    v = v0;
    tf = Module::getImplementation<IFreezingPoint>()(sss[0]);
    cpml = Water::rho * Water::cp * mld[0];
    qio = qio0;

    /* It's only the SSH gradient which has an effect, so being able to sett a constant SSH is
     * useless. */
    ssh = 0.;
}

UniformOcean& UniformOcean::setSST(double sstIn)
{
    sst0 = sstIn;
    return *this;
}
UniformOcean& UniformOcean::setSSS(double sssIn)
{
    sss0 = sssIn;
    return *this;
}
UniformOcean& UniformOcean::setMLD(double mldIn)
{
    mld0 = mldIn;
    return *this;
}
UniformOcean& UniformOcean::setU(double uIn)
{
    u0 = uIn;
    return *this;
}
UniformOcean& UniformOcean::setV(double vIn)
{
    v0 = vIn;
    return *this;
}
UniformOcean& UniformOcean::setQio(double qioIn)
{
    qio0 = qioIn;
    return *this;
}

} /* namespace Nextsim */
