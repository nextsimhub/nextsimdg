/*!
 * @file IIceThermodynamics.hpp
 *
 * @date Mar 16, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IICETHERMODYNAMICS_HPP
#define IICETHERMODYNAMICS_HPP

#include "include/ModelArray.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/Time.hpp"

namespace Nextsim {
//! An interface class to update the ice thermodynamics.
class IIceThermodynamics : public ModelComponent {
public:
    ~IIceThermodynamics() = default;

    std::string getName() const override { return "IceThermodynamics"; }
    void setData(const ModelState::DataMap& ms) override
    {
        tice.resize();
        deltaHi.resize();
        snowToIce.resize();
    }
    ModelState getState() const override { return ModelState(); }
    ModelState getState(const OutputLevel&) const override { return getState(); }
    ModelState getStateRecursive(const OutputSpec& os) const override
    {
        return os ? getState() : ModelState();
    }
    /*!
     * Updates the ice thermodynamic and thickness growth calculation for the timestep.
     *
     * @param tStep The object containing the timestep start and duration times.
     */
    virtual void update(const TimestepTime& tsTime) = 0;

protected:
    IIceThermodynamics()
        : tice(ModelArray::Type::Z)
        , deltaHi(ModelArray::Type::H)
        , snowToIce(ModelArray::Type::H)
    {
        registerModule();

        ModelComponent::registerSharedArray(SharedArray::DELTA_HICE, &deltaHi);
        ModelComponent::registerSharedArray(SharedArray::T_ICE, &tice);
    }

    ModelArrayRef<SharedArray::H_ICE, RW> hice; // From IceGrowth
    ModelArrayRef<SharedArray::C_ICE, RW> cice; // From IceGrowth
    ModelArrayRef<SharedArray::H_SNOW, RW> hsnow; // From Ice Growth
    ModelArrayRef<SharedArray::Q_IC, RW>
        qic; // From IceTemperature. Conductive heat flux to the ice surface.
    ModelArrayRef<SharedArray::Q_IO, RW> qio; // From FluxCalculation
    ModelArrayRef<SharedArray::Q_IA, RO> qia; // From FluxCalculation
    ModelArrayRef<SharedArray::DQIA_DT, RO> dQia_dt; // From FluxCalculation
    ModelArrayRef<SharedArray::SUBLIM, RO> sublim; // From AtmosphereState
    ModelArrayRef<ProtectedArray::T_ICE> tice0; // Timestep initial ice temperature
    ModelArrayRef<ProtectedArray::TF> tf; // Sea water freezing temperature
    ModelArrayRef<ProtectedArray::SNOW> snowfall; // From ExternalData
    ModelArrayRef<ProtectedArray::SSS> sss; // From ExternalData (possibly PrognosticData)
    // Owned, shared arrays
    HField tice;
    HField deltaHi;
    // Owned, Module-private arrays
    HField snowToIce;
};

} /* namespace Nextsim */

#endif /* IICETHERMODYNAMICS_HPP */
