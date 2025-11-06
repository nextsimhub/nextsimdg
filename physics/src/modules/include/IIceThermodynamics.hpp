/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IICETHERMODYNAMICS_HPP
#define IICETHERMODYNAMICS_HPP

#include "include/ConfigurationHelp.hpp"
#include "include/FieldAdvection.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelArrayAccessor.hpp"
#include "include/ModelArraySlice.hpp"
#include "include/ModelComponent.hpp"
#include "include/Slice.hpp"
#include "include/Time.hpp"
#include "include/gridNames.hpp"

namespace Nextsim {
//! An interface class to update the ice thermodynamics.
class IIceThermodynamics : public ModelComponent {
public:
    ~IIceThermodynamics() = default;

    std::string getName() const override { return "IceThermodynamics"; }
    void setData(const ModelState::DataMap& ms) override
    {
        deltaHiAccessor.getDeviceRO();
        tsurf.resize();
        deltaHi.resize();
        snowToIce.resize();

        /* If the surface temperature is not in the restart file, then we simply set it to the
         * zero. It's a safe approximation, and it seems the user doesn't
         * really care! */
        try {
            tsurf = ms.at(tsurfName);
        } catch (const std::out_of_range& e) {
            Logged::info("No " + tsurfName + " field in restart file. Setting it to 0°C.\n");
            tsurf = 0.;
        }
    }

    ModelState getStatePrognostic() const override
    {
        return { {
                     { tsurfName, tsurf },
                 },
            getConfiguration() };
    }

    ModelState getStateDiagnostic() const override
    {
        ModelState state = { {
                                 { "delta_H_ice", deltaHi },
                                 { "snow_to_ice", snowToIce },
                             },
            getConfiguration() };
        state.merge(getStatePrognostic());

        return state;
    }

    /*!
     * Updates the ice thermodynamic and thickness growth calculation for the timestep.
     *
     * @param tStep The object containing the timestep start and duration times.
     */
    virtual void update(const TimestepTime& tsTime)
    {
        FieldAdvection::advectField(tsurf, tsTime, minT, 0.);
    }

    inline static std::string getKappaSConfigKey() { return "nextsim_thermo.ks"; }

protected:
    IIceThermodynamics()
        : tsurf(ModelArray::AdvectionType)
        , deltaHi(ModelArray::Type::H)
        , snowToIce(ModelArray::Type::H)
        , hice(getStore())
        , cice(getStore())
        , hsnow(getStore())
        , qic(getStore())
        , qio(getStore())
        , qow(getStore())
        , qia(getStore())
        , dQia_dt(getStore())
        , penSw(getStore())
        , sublim(getStore())
        , tf(getStore())
        , snowfall(getStore())
        , sss(getStore())
        , qswBase(getStore())
        // proposed interface
        , hiceAccessor(getStore())
        , qiaAccessor(getStore())
        // formerly owned arrays are initialized by special constructor
        , deltaHiAccessor(getStore(), RW, ModelArray::Type::H)
        //, tsurfAccessor(getStore(), RO, ModelArray::Type::H)
    {
        getStore().registerArray(Shared::DELTA_HICE, &deltaHi, RW);
        getStore().registerArray(Protected::T_SURF, &tsurf, RO);
    }

    ModelArrayRef<Shared::H_ICE_DG, RW> hice; // From PrognosticData
    ModelArrayRef<Shared::C_ICE_DG, RW> cice; // From PrognosticData
    ModelArrayRef<Shared::H_SNOW_DG, RW> hsnow; // From PrognosticData
    ModelArrayRef<Shared::Q_IC, RW>
        qic; // From IceTemperature. Conductive heat flux to the ice surface.
    ModelArrayRef<Shared::Q_SW_BASE, RW> qswBase; // Short-wave flux through the base of the ice
    ModelArrayRef<Shared::Q_IO, RW> qio; // From FluxCalculation
    ModelArrayRef<Shared::Q_OW, RW> qow; // From FluxCalculation
    ModelArrayRef<Shared::Q_IA, RO> qia; // From FluxCalculation
    ModelArrayRef<Shared::DQIA_DT, RO> dQia_dt; // From FluxCalculation
    ModelArrayRef<Shared::Q_PEN_SW, RO> penSw; // From FluxCalculation
    ModelArrayRef<Shared::SUBLIM, RO> sublim; // From AtmosphereState
    ModelArrayRef<Protected::TF> tf; // Sea water freezing temperature
    ModelArrayRef<Protected::SNOW> snowfall; // From ExternalData
    ModelArrayRef<Protected::SSS> sss; // From ExternalData (possibly PrognosticData)
    // Owned, shared arrays
    AdvectedField tsurf;
    HField deltaHi;
    // Owned, Module-private arrays
    HField snowToIce;

    // proposed interface
    ModelArrayAccessor<Shared::H_ICE_DG, RW> hiceAccessor;
    ModelArrayAccessor<Shared::Q_IA, RO> qiaAccessor;
    ModelArrayAccessor<Shared::DELTA_HICE, RW> deltaHiAccessor;
    //ModelArrayAccessor<Shared::T_SURF, RW> tsurfAccessor;

    constexpr static double minT = -90.0;
};

} /* namespace Nextsim */

#endif /* IICETHERMODYNAMICS_HPP */
