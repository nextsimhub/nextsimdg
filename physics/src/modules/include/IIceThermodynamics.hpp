/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IICETHERMODYNAMICS_HPP
#define IICETHERMODYNAMICS_HPP

#include "include/ConfigurationHelp.hpp"
#include "include/FieldAdvection.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArrayAccessor.hpp"
#include "include/ModelArraySlice.hpp"
#include "include/ModelComponent.hpp"
#include "include/Slice.hpp"
#include "include/Time.hpp"
#ifdef USE_XIOS
#include "include/Xios.hpp"
#endif
#include "include/gridNames.hpp"

namespace Nextsim {

namespace Private {
    inline constexpr TextTag SNOW_TO_ICE = "SNOW_TO_ICE";
}

//! An interface class to update the ice thermodynamics.
class IIceThermodynamics : public ModelComponent {
public:
    ~IIceThermodynamics() = default;

    std::string getName() const override { return "IceThermodynamics"; }
    void setData(const ModelState::DataMap& ms) override
    {
        AdvectedField& tsurf = tsurfAccessor.getHostRW();
        tsurf.reinitialize();
        HField& deltaHi = deltaHiAccessor.getHostRW();
        deltaHi.reinitialize();
        snowToIceAccessor.getHostRW().reinitialize();

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
                     { tsurfName, tsurfAccessor.getHostRO() },
                 },
            getConfiguration() };
    }

    ModelState getStateDiagnostic() const override
    {
        ModelState state = { {
                                 { "delta_H_ice", deltaHiAccessor.getHostRO() },
                                 { "snow_to_ice", snowToIceAccessor.getHostRO() },
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
#ifdef USE_KOKKOS
        FieldAdvection::advectField(tsurfAccessor.getDeviceRW(), tsTime, minT, FloatType(0));
#else
        FieldAdvection::advectField(tsurfAccessor.getHostRW(), tsTime, minT, FloatType(0));
#endif
    }

    inline static std::string getKappaSConfigKey() { return "nextsim_thermo.ks"; }

protected:
    IIceThermodynamics()
        : tsurfAccessor(getStore(), RO, ModelArray::AdvectionType)
        , deltaHiAccessor(getStore(), RW, ModelArray::Type::H)
        , snowToIceAccessor(getStore(), RO, ModelArray::Type::H)
        , hiceAccessor(getStore())
        , ciceAccessor(getStore())
        , hsnowAccessor(getStore())
        //    , qicAccessor(getStore())
        , qioAccessor(getStore())
        , qowAccessor(getStore())
        , qiaAccessor(getStore())
        , dQia_dtAccessor(getStore())
        , penSwAccessor(getStore())
        , sublimAccessor(getStore())
        , tfAccessor(getStore())
        , snowfallAccessor(getStore())
        , sssAccessor(getStore())
        , qswBaseAccessor(getStore())
    {

#ifdef USE_XIOS
        // Set XIOS field types
        Xios& xiosHandler = Xios::getInstance();
        xiosHandler.setPrognosticFieldType(tsurfName, ModelArray::AdvectionType);
#endif
    }

    ModelArrayAccessor<Shared::H_ICE_DG, RW> hiceAccessor; // From PrognosticData
    ModelArrayAccessor<Shared::C_ICE_DG, RW> ciceAccessor; // From PrognosticData
    ModelArrayAccessor<Shared::H_SNOW_DG, RW> hsnowAccessor; // From PrognosticData
    // Q_IC is owned by ThermoIce0 and currently not used
    // ModelArrayAccessor<Shared::Q_IC, RW>
    //    qicAccessor; // From IceTemperature. Conductive heat flux to the ice surface.
    ModelArrayAccessor<Shared::Q_SW_BASE, RW>
        qswBaseAccessor; // Short-wave flux through the base of the ice
    ModelArrayAccessor<Shared::Q_IO, RW> qioAccessor; // From FluxCalculation
    ModelArrayAccessor<Shared::Q_OW, RW> qowAccessor; // From FluxCalculation
    ModelArrayAccessor<Shared::Q_IA, RO> qiaAccessor; // From FluxCalculation
    ModelArrayAccessor<Shared::DQIA_DT, RO> dQia_dtAccessor; // From FluxCalculation
    ModelArrayAccessor<Shared::Q_PEN_SW, RO> penSwAccessor; // From FluxCalculation
    ModelArrayAccessor<Shared::SUBLIM, RO> sublimAccessor; // From AtmosphereState
    ModelArrayAccessor<Protected::TF> tfAccessor; // Sea water freezing temperature
    ModelArrayAccessor<Protected::SNOW> snowfallAccessor; // From ExternalData
    ModelArrayAccessor<Protected::SSS> sssAccessor; // From ExternalData (possibly PrognosticData)
    // Owned, shared arrays
    ModelArrayAccessor<Protected::T_SURF, RW> tsurfAccessor;
    ModelArrayAccessor<Shared::DELTA_HICE, RW> deltaHiAccessor;
    // Owned, Module-private arrays
    ModelArrayAccessor<Private::SNOW_TO_ICE, RW> snowToIceAccessor;

    constexpr static FloatType minT = -90.0;
};

} /* namespace Nextsim */

#endif /* IICETHERMODYNAMICS_HPP */
