/*!
 * @file IIceOceanHeatFluxModule.cpp
 *
 * @date Feb 21, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include <string>
#include "include/BasicIceOceanHeatFlux.hpp"
#include "include/IceOceanHeatFluxModule.hpp"

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
    = std::move(newImpl<Nextsim::IIceOceanHeatFlux, Nextsim::BasicIceOceanHeatFlux>());

template <> std::string Module<Nextsim::IIceOceanHeatFlux>::moduleName()
{
    return "IIceOceanHeatFlux";
}

template <> Nextsim::IIceOceanHeatFlux& getImplementation<Nextsim::IIceOceanHeatFlux>()
{
    return getImplTemplate<Nextsim::IIceOceanHeatFlux, IceOceanHeatFluxModule>();
}
template <> void setImplementation<Nextsim::IIceOceanHeatFlux>(const std::string& implName)
{
    setImplTemplate<IceOceanHeatFluxModule>(implName);
}
template <> std::unique_ptr<Nextsim::IIceOceanHeatFlux> getInstance()
{
    return getInstTemplate<Nextsim::IIceOceanHeatFlux, IceOceanHeatFluxModule>();
}
} /* namespace Module */
