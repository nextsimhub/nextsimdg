/*!
 * @file IVerticalIceGrowth.hpp
 *
 * @date Mar 16, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IVERTICALICEGROWTH_HPP
#define IVERTICALICEGROWTH_HPP

#include "include/ModelArray.hpp"
#include "include/ModelModule.hpp"

namespace Nextsim {
class IVerticalIceGrowth : public ModelModule {
public:
    ~IVerticalIceGrowth() = default;

protected:
    IVerticalIceGrowth()
    {
        ModelModule::requestSharedArray(SharedArray::H_ICE, &hice);
        ModelModule::requestSharedArray(SharedArray::C_ICE, &cice);
        ModelModule::requestSharedArray(SharedArray::H_SNOW, &hsnow);
        ModelModule::requestSharedArray(SharedArray::T_ICE, &tice);
        ModelModule::requestSharedArray(SharedArray::Q_IA, &qia);
        ModelModule::requestSharedArray(SharedArray::Q_IC, &qic);
        ModelModule::requestSharedArray(SharedArray::Q_IO, &qio);
        ModelModule::requestSharedArray(SharedArray::DQIA_DT, &dQia_dT);
        ModelModule::requestSharedArray(SharedArray::SUBLIM, &sublim);
        ModelModule::requestProtectedArray(ProtectedArray::SNOW, &snowfall);
    }

    HField* hice; // From IceGrowth
    HField* cice; // From IceGrowth
    HField* hsnow; // From Ice Growth?
    ZField* tice; // Updated value from IceTemperature
    HField* qic; // From IceTemperature. Conductive heat flux to the ice surface.
    HField* qia; // From FluxCalculation
    HField* dQia_dT; // From FluxCalculation
    HField* qio; // From FluxCalculation
    HField* sublim; // From AtmosphereState
    const HField* snowfall; // From ExternalData
    // And owns nothing of its own
};

} /* namespace Nextsim */

#endif /* IVERTICALICEGROWTH_HPP */
