/*!
 * @file IIceThermodynamics.hpp
 *
 * @date Mar 16, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IICETHERMODYNAMICS_HPP
#define IICETHERMODYNAMICS_HPP

#include "include/ModelArrayRef.hpp"
#include "include/ConfigurationHelp.hpp"
#include "include/ModelArray.hpp"
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

    virtual void initialiseTice() { tice = tice0; }

    inline static std::string getKappaSConfigKey() { return "nextsim_thermo.ks"; }

    virtual size_t getNZLevels() const = 0;

protected:
    IIceThermodynamics()
        : tice(ModelArray::Type::Z)
        , deltaHi(ModelArray::Type::H)
        , snowToIce(ModelArray::Type::H)
        , hice(getStore())
        , cice(getStore())
        , hsnow(getStore())
        , qic(getStore())
        , qio(getStore())
        , qia(getStore())
        , dQia_dt(getStore())
        , penSw(getStore())
        , sublim(getStore())
        , tice0(getStore())
        , tf(getStore())
        , snowfall(getStore())
        , sss(getStore())
    {
        getStore().registerArray(Shared::DELTA_HICE, &deltaHi, RW);
        getStore().registerArray(Shared::T_ICE, &tice, RW);
    }

    ModelArrayRef<Shared::H_ICE, RW> hice; // From IceGrowth
    ModelArrayRef<Shared::C_ICE, RW> cice; // From IceGrowth
    ModelArrayRef<Shared::H_SNOW, RW> hsnow; // From Ice Growth
    ModelArrayRef<Shared::Q_IC, RW>
        qic; // From IceTemperature. Conductive heat flux to the ice surface.
    ModelArrayRef<Shared::Q_IO, RW> qio; // From FluxCalculation
    ModelArrayRef<Shared::Q_IA, RO> qia; // From FluxCalculation
    ModelArrayRef<Shared::DQIA_DT, RO> dQia_dt; // From FluxCalculation
    ModelArrayRef<Shared::Q_PEN_SW, RO> penSw; // From FluxCalculation
    ModelArrayRef<Shared::SUBLIM, RO> sublim; // From AtmosphereState
    ModelArrayRef<Protected::T_ICE> tice0; // Timestep initial ice temperature
    ModelArrayRef<Protected::TF> tf; // Sea water freezing temperature
    ModelArrayRef<Protected::SNOW> snowfall; // From ExternalData
    ModelArrayRef<Protected::SSS> sss; // From ExternalData (possibly PrognosticData)
    // Owned, shared arrays
    HField tice;
    HField deltaHi;
    // Owned, Module-private arrays
    HField snowToIce;
};

} /* namespace Nextsim */

#endif /* IICETHERMODYNAMICS_HPP */
