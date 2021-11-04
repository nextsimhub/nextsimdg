/*!
 * @file Physics1dBase.cpp
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/Physics1dBase.hpp"

namespace Nextsim {

Physics1dBase::Physics1dBase()
{
    // TODO Auto-generated constructor stub
}

Physics1dBase::~Physics1dBase()
{
    // TODO Auto-generated destructor stub
}

template <class Phys> void Physics1dBase::physics1d(ElementData<Phys>& data)
{
    Phys::updateDerivedData(data);
    Phys::massFluxOpenWater(data);
    Phys::momentumFluxOpenWater(data);
    Phys::heatFluxOpenWater(data);

    Phys::massFluxIceAtmosphere(data);
    // Ice momentum fluxes are handled by the dynamics
    Phys::heatFluxIceAtmosphere(data);

    Phys::massFluxIceOcean(data);
    // Ice momentum fluxes are handled by the dynamics
    Phys::heatFluxIceOcean(data);
}
} /* namespace Nextsim */
