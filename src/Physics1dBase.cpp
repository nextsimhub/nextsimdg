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

template <class T>
void Physics1dBase::physics1d<T>(ElementData& data) {
    T::updateDerivedData(data);
    T::massFluxOpenWater(data);
    T::momentumFluxOpenWater(data);
    T::heatFluxOpenWater(data);

    T::massFluxIceAtmosphere(data);
    // Ice momentum fluxes are handled by the dynamics
    T::heatFluxIceAtmosphere(data);

    T::massFluxIceOcean(data);
    // Ice momentum fluxes are handled by the dynamics
    T::heatFluxIceOcean(data);
}
} /* namespace Nextsim */
