/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/UniformOcean.hpp"

#include "include/IFreezingPoint.hpp"
#include "include/NextsimModule.hpp"
#include "include/constants.hpp"

namespace Nextsim {

void UniformOcean::setData(const ModelState::DataMap& ms)
{
    IOceanBoundary::setData(ms);
    sstAccessor.getHostRW() = sst0;
    HField& sss = sssAccessor.getHostRW();
    sss = sss0;
    HField& mld = mldAccessor.getHostRW();
    mld = mld0;
    uAccessor.getHostRW() = u0;
    vAccessor.getHostRW() = v0;
    tfAccessor.getHostRW() = Module::getImplementation<IFreezingPoint>()(sss[0]);
    cpmlAccessor.getHostRW() = Water::rho * Water::cp * mld[0];
    qioAccessor.getHostRW() = qio0;

    /* It's only the SSH gradient which has an effect, so being able to set a constant SSH is
     * useless. */
    sshAccessor.getHostRW() = 0.;
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
