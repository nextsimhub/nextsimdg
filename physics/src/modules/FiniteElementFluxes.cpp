/*!
 * @file FiniteElementFluxes.cpp
 *
 * @date Apr 29, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/FiniteElementFluxes.hpp"

namespace Nextsim {

ModelState FiniteElementFluxes::getState() const
{
    return ModelState();
}
ModelState FiniteElementFluxes::getState(const OutputLevel&) const
{
    return getState();
}

void FiniteElementFluxes::update(const TimestepTime& tst)
{

}

} /* namespace Nextsim */
