/*!
 * @file NoCoupler.cpp
 *
 * @date Jul 22, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/Coupler.hpp"

namespace Nextsim {

//! The Coupler class implementation when there is no coupler: configuration.
void Coupler::configure() { }
//! The Coupler class implementation when there is no coupler: initialization.
void Coupler::initialize(const ModelMetadata& meta) { }
//! The Coupler class implementation when there is no coupler: shut down.
void Coupler::terminate() { }

} /* namespace Nextsim */
