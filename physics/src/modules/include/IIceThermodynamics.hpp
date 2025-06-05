/*!
 * @file IIceThermodynamics.hpp
 *
 * @date 24 Sep 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IICETHERMODYNAMICS_HPP
#define IICETHERMODYNAMICS_HPP

#include "include/ConfigurationHelp.hpp"
#include "include/gridNames.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelArraySlice.hpp"
#include "include/ModelComponent.hpp"
#include "include/Slice.hpp"
#include "include/Time.hpp"

namespace Nextsim {
//! An interface class to update the ice thermodynamics.
class IIceThermodynamics : public ModelComponent {
public:
    ~IIceThermodynamics() = default;

    std::string getName() const override { return "IceThermodynamics"; }
    void setData(const ModelState::DataMap& ms) override
    {
        tsurf.resize();
        if (ms.count(tsurfName) > 0) {
            tsurf = ms.at(tsurfName);
        } else if (static_cast<ModelArray>(tice0).nDimensions() == 2) {
            tsurf = tice0;
        } else {
            tsurf = static_cast<ModelArray>(tice0)[z0Slice];
        }
        deltaHi.resize();
        snowToIce.resize();
    }

    ModelState getStatePrognostic() const override {
        return { {
            { tsurfName, tsurf },
        }, getConfiguration() };
    }

    ModelState getStateDiagnostic() const override
    {
        ModelState state = { {
            { "delta_H_ice", deltaHi },
            { "snow_to_ice", snowToIce },
        },
                getConfiguration()
        };
        state.merge(getStatePrognostic());

        return state;
    }

    /*!
     * Updates the ice thermodynamic and thickness growth calculation for the timestep.
     *
     * @param tStep The object containing the timestep start and duration times.
     */
    virtual void update(const TimestepTime& tsTime) = 0;

    inline static std::string getKappaSConfigKey() { return "nextsim_thermo.ks"; }

protected:
    IIceThermodynamics()
        : tsurf(ModelArray::Type::H)
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
        , qswBase(getStore())
    {
        getStore().registerArray(Shared::DELTA_HICE, &deltaHi, RW);
        getStore().registerArray(Protected::T_SURF, &tsurf, RO);
    }

    ModelArrayRef<Shared::H_ICE, RW> hice; // From IceGrowth
    ModelArrayRef<Shared::C_ICE, RW> cice; // From IceGrowth
    ModelArrayRef<Shared::H_SNOW, RW> hsnow; // From Ice Growth
    ModelArrayRef<Shared::Q_IC, RW>
        qic; // From IceTemperature. Conductive heat flux to the ice surface.
    ModelArrayRef<Shared::Q_SW_BASE, RW> qswBase; // Short-wave flux through the base of the ice
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
    HField tsurf;
    HField deltaHi;
    // Owned, Module-private arrays
    HField snowToIce;

    const ArraySlicer::Slice z0Slice {{{ }, { }, {0}}};
};

} /* namespace Nextsim */

#endif /* IICETHERMODYNAMICS_HPP */
