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
#include "include/Time.hpp"

namespace Nextsim {
class IVerticalIceGrowth : public ModelModule {
public:
    ~IVerticalIceGrowth() = default;

    std::string getName() const override { return "VerticalIceGrowth"; }
    void setData(const ModelState& ms) override { }
    ModelState getState() const override { return ModelState(); }
    ModelState getState(const OutputLevel&) const override { return getState(); }
    virtual void update(const TimestepTime& tsTime) = 0;

protected:
    IVerticalIceGrowth()
    {
        registerModule();
        ModelModule::requestSharedArray(SharedArray::H_ICE, &hice);
        ModelModule::requestSharedArray(SharedArray::C_ICE, &cice);
        ModelModule::requestSharedArray(SharedArray::H_SNOW, &hsnow);
        ModelModule::requestSharedArray(SharedArray::T_ICE, &tice);
        ModelModule::requestSharedArray(SharedArray::Q_IC, &qic);
        ModelModule::requestSharedArray(SharedArray::Q_IO, &qio);
        ModelModule::requestProtectedArray(SharedArray::Q_IA, &qia);
        ModelModule::requestProtectedArray(SharedArray::DQIA_DT, &dQia_dT);
        ModelModule::requestProtectedArray(SharedArray::SUBLIM, &sublim);
        ModelModule::requestProtectedArray(ProtectedArray::SNOW, &snowfall);
    }

    HField* hice; // From IceGrowth
    HField* cice; // From IceGrowth
    HField* hsnow; // From Ice Growth?
    ZField* tice; // Updated value from IceTemperature
    HField* qic; // From IceTemperature. Conductive heat flux to the ice surface.
    HField* qio; // From FluxCalculation
    const HField* qia; // From FluxCalculation
    const HField* dQia_dT; // From FluxCalculation
    const HField* sublim; // From AtmosphereState
    const HField* snowfall; // From ExternalData
    // Owned, Module-private arrays
    HField snowToIce;
};

} /* namespace Nextsim */

#endif /* IVERTICALICEGROWTH_HPP */
