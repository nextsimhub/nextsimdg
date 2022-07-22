/*!
 * @file OasisCoupler.cpp
 *
 * @date Jul 22, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/Coupler.hpp"
#include "include/OasisData.hpp"
#include "oasis_c.h"

namespace Nextsim {
//! The OASIS Coupler class implementation: configuration.
void Coupler::configure()
{
    pData = new CouplerData;
}
//! The OASIS Coupler class implementation: initialization.
void Coupler::initialize(const ModelMetadata& meta)
{
    // Will initialize the Oasis data class for this process

}
//! The OASIS Coupler class implementation: shut down.
void Coupler::terminate()
{
    delete pData;
}

} /* namespace Nextsim */
