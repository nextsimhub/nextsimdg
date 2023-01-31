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
        , hice(getSharedArray())
        , cice(getSharedArray())
        , hsnow(getSharedArray())
        , qic(getSharedArray())
        , qio(getSharedArray())
        , qia(getSharedArray())
        , dQia_dt(getSharedArray())
        , sublim(getSharedArray())
        , tice0(getProtectedArray())
        , tf(getProtectedArray())
        , snowfall(getProtectedArray())
        , sss(getProtectedArray())
    {
        registerModule();

        ModelComponent::registerSharedArray(SharedArray::DELTA_HICE, &deltaHi);
        ModelComponent::registerSharedArray(SharedArray::T_ICE, &tice);
    }

    ModelArrayRef<SharedArray::H_ICE, MARBackingStore, RW> hice; // From IceGrowth
    ModelArrayRef<SharedArray::C_ICE, MARBackingStore, RW> cice; // From IceGrowth
    ModelArrayRef<SharedArray::H_SNOW, MARBackingStore, RW> hsnow; // From Ice Growth
    ModelArrayRef<SharedArray::Q_IC, MARBackingStore, RW>
        qic; // From IceTemperature. Conductive heat flux to the ice surface.
    ModelArrayRef<SharedArray::Q_IO, MARBackingStore, RW> qio; // From FluxCalculation
    ModelArrayRef<SharedArray::Q_IA, MARBackingStore, RO> qia; // From FluxCalculation
    ModelArrayRef<SharedArray::DQIA_DT, MARBackingStore, RO> dQia_dt; // From FluxCalculation
    ModelArrayRef<SharedArray::SUBLIM, MARBackingStore, RO> sublim; // From AtmosphereState
    ModelArrayRef<ProtectedArray::T_ICE, MARConstBackingStore>
        tice0; // Timestep initial ice temperature
    ModelArrayRef<ProtectedArray::TF, MARConstBackingStore> tf; // Sea water freezing temperature
    ModelArrayRef<ProtectedArray::SNOW, MARConstBackingStore> snowfall; // From ExternalData
    ModelArrayRef<ProtectedArray::EXT_SSS, MARConstBackingStore>
        sss; // From ExternalData (possibly PrognosticData)
    // Owned, shared arrays
    HField tice;
    HField deltaHi;
    // Owned, Module-private arrays
    HField snowToIce;
};

} /* namespace Nextsim */

#endif /* IICETHERMODYNAMICS_HPP */
