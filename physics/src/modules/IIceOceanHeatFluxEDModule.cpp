/*!
 * @file IIceOceanHeatFluxModule.cpp
 *
 * @date Feb 21, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include <string>
#include "include/BasicIceOceanHeatFluxED.hpp"
#include "include/IIceOceanHeatFluxEDModule.hpp"

namespace Module {
const std::string BASICICEOCEANHEATFLUX = "BasicIceOceanHeatFlux";

template <>
Module<Nextsim::IIceOceanHeatFluxED>::map Module<Nextsim::IIceOceanHeatFluxED>::functionMap = {
    { BASICICEOCEANHEATFLUX, newImpl<Nextsim::IIceOceanHeatFluxED, Nextsim::BasicIceOceanHeatFluxED> },
};

template <>
Module<Nextsim::IIceOceanHeatFluxED>::fn Module<Nextsim::IIceOceanHeatFluxED>::spf
    = functionMap.at(BASICICEOCEANHEATFLUX);
template <>
std::unique_ptr<Nextsim::IIceOceanHeatFluxED> Module<Nextsim::IIceOceanHeatFluxED>::staticInstance
    = std::move(newImpl<Nextsim::IIceOceanHeatFluxED, Nextsim::BasicIceOceanHeatFluxED>());

template <> std::string Module<Nextsim::IIceOceanHeatFluxED>::moduleName()
{
    return "IIceOceanHeatFlux";
}

template <> Nextsim::IIceOceanHeatFluxED& getImplementation<Nextsim::IIceOceanHeatFluxED>()
{
    return getImplTemplate<Nextsim::IIceOceanHeatFluxED, IIceOceanHeatFluxEDModule>();
}
template <> void setImplementation<Nextsim::IIceOceanHeatFluxED>(const std::string& implName)
{
    setImplTemplate<IIceOceanHeatFluxEDModule>(implName);
}
template <> std::unique_ptr<Nextsim::IIceOceanHeatFluxED> getInstance()
{
    return getInstTemplate<Nextsim::IIceOceanHeatFluxED, IIceOceanHeatFluxEDModule>();
}
} /* namespace Module */
