/*!
 * @file IIceOceanHeatFluxModule.cpp
 *
 * @date Feb 21, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/IIceOceanHeatFluxModule.hpp"

#include "include/BasicIceOceanHeatFlux.hpp"

#include <string>

namespace Module {
const std::string BASICICEOCEANHEATFLUX = "BasicIceOceanHeatFlux";

template <>
Module<Nextsim::IIceOceanHeatFlux>::map Module<Nextsim::IIceOceanHeatFlux>::functionMap = {
    { BASICICEOCEANHEATFLUX, newImpl<Nextsim::IIceOceanHeatFlux, Nextsim::BasicIceOceanHeatFlux> },
};

template <>
Module<Nextsim::IIceOceanHeatFlux>::fn Module<Nextsim::IIceOceanHeatFlux>::spf
    = functionMap.at(BASICICEOCEANHEATFLUX);
template <>
std::unique_ptr<Nextsim::IIceOceanHeatFlux> Module<Nextsim::IIceOceanHeatFlux>::staticInstance
    = std::move(Module<Nextsim::IIceOceanHeatFlux>::spf());

template <> std::string Module<Nextsim::IIceOceanHeatFlux>::moduleName()
{
    return "IIceOceanHeatFlux";
}

template <> Nextsim::IIceOceanHeatFlux& getImplementation<Nextsim::IIceOceanHeatFlux>()
{
    return getImplTemplate<Nextsim::IIceOceanHeatFlux, IIceOceanHeatFluxModule>();
}
template <> void setImplementation<Nextsim::IIceOceanHeatFlux>(const std::string& implName)
{
    setImplTemplate<IIceOceanHeatFluxModule>(implName);
}
template <> std::unique_ptr<Nextsim::IIceOceanHeatFlux> getInstance()
{
    return getInstTemplate<Nextsim::IIceOceanHeatFlux, IIceOceanHeatFluxModule>();
}
} /* namespace Module */
