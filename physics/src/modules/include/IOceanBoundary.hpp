/*!
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef IOCEANBOUNDARY_HPP
#define IOCEANBOUNDARY_HPP

#include "include/CheckingModelComponent.hpp"
#include "include/KernelAlternatives.hpp"
#include "include/ModelArrayAccessor.hpp"
#include "include/constants.hpp"
#include "include/gridNames.hpp"

namespace Nextsim {

//! An interface class for the oceanic inputs into the ice physics.
class IOceanBoundary : public CheckingModelComponent {
public:
    IOceanBoundary()
        : qioAccessor(getStore(), RW, ModelArray::Type::H)
        , sstAccessor(getStore(), RO, ModelArray::Type::H, std::pair(-5.0, 50.0))
        , sssAccessor(getStore(), RO, ModelArray::Type::H, std::pair(0.0, 50.0))
        , mldAccessor(getStore(), RO, ModelArray::Type::H, std::pair(1e-3, 12e3))
        , cpmlAccessor(getStore(), RO, ModelArray::Type::H, std::pair(0.0, 1e11))
        , tfAccessor(getStore(), RO, ModelArray::Type::H, std::pair(-5.0, 0.0))
        , uAccessor(getStore(), RO, ModelArray::Type::H, std::pair(-1.0, 1.0))
        , vAccessor(getStore(), RO, ModelArray::Type::H, std::pair(-1.0, 1.0))
        , sshAccessor(getStore(), RO, ModelArray::Type::H, std::pair(-10.0, 10.0))
        , qNoSunAccessor(m_couplingArrays, RO, ModelArray::Type::H)
        , qswNetAccessor(m_couplingArrays, RO, ModelArray::Type::H, std::pair(-1e3, 1e3))
        , fwFluxAccessor(m_couplingArrays, RO, ModelArray::Type::H)
        , sFluxAccessor(m_couplingArrays, RO, ModelArray::Type::H)
        , qswowAccessor(getStore(), RW, ModelArray::Type::H, std::pair(-1e3, 1e-3))
        , qswBaseAccessor(getStore(), RW, ModelArray::Type::H, std::pair(-1e3, 1e-6))
        , tauXAccessor(m_couplingArrays, RO, ModelArray::Type::H, std::pair(-10.0, 10.0))
        , tauYAccessor(m_couplingArrays, RO, ModelArray::Type::H, std::pair(-10.0, 10.0))
        , ciceAccessor(getStore())
        , evapAccessor(getStore())
        , rainAccessor(getStore())
        , newIceAccessor(getStore())
        , deltaHiceAccessor(getStore())
        , deltaSmeltAccessor(getStore())
        , qowAccessor(getStore())
        , tauXIOAccessor(getStore())
        , tauYIOAccessor(getStore())
        , tauXOWAccessor(getStore())
        , tauYOWAccessor(getStore())
    {
    }
    virtual ~IOceanBoundary() = default;

    std::string getName() const override { return "IOceanBoundary"; }
    void setData(const ModelState::DataMap& ms) override
    {
        HField& qio = qioAccessor.getHostRW();
        qio.reinitialize();
        HField& sst = sstAccessor.getHostRW();
        sst.reinitialize();
        HField& sss = sssAccessor.getHostRW();
        sss.reinitialize();
        HField& mld = mldAccessor.getHostRW();
        mld.reinitialize();
        HField& cpml = cpmlAccessor.getHostRW();
        cpml.reinitialize();
        HField& tf = tfAccessor.getHostRW();
        tf.reinitialize();
        UField& u = uAccessor.getHostRW();
        u.reinitialize();
        VField& v = vAccessor.getHostRW();
        v.reinitialize();
        HField& ssh = sshAccessor.getHostRW();
        ssh.reinitialize();
        HField& qNoSun = qNoSunAccessor.getHostRW();
        qNoSun.reinitialize();
        HField& qswNet = qswNetAccessor.getHostRW();
        qswNet.reinitialize();
        HField& fwFlux = fwFluxAccessor.getHostRW();
        fwFlux.reinitialize();
        HField& sFlux = sFluxAccessor.getHostRW();
        sFlux.reinitialize();
        HField& qswow = qswowAccessor.getHostRW();
        qswow.reinitialize();
        HField& qswBase = qswBaseAccessor.getHostRW();
        qswBase.reinitialize();
        HField& tauX = tauXAccessor.getHostRW();
        tauX.reinitialize();
        HField& tauY = tauYAccessor.getHostRW();
        tauY.reinitialize();

        addChecks({
            { "qio", qioAccessor },
            { "sst", sstAccessor },
            { "sss", sssAccessor },
            { "mld", mldAccessor },
            { "cpml", cpmlAccessor },
            { "tf", tfAccessor },
            { "u", uAccessor },
            { "v", vAccessor },
            { "ssh", sshAccessor },
            { "qNoSun", qNoSunAccessor },
            { "qswNet", qswNetAccessor },
            { "fwFlux", fwFluxAccessor },
            { "sFlux", sFluxAccessor },
            { "qswow", qswowAccessor },
            { "qswBase", qswBaseAccessor },
            { "tauX", tauXAccessor },
            { "tauY", tauYAccessor },
        });

        if (ms.count(sstName)) {
            sst = ms.at(sstName);
        }
        if (ms.count(sssName)) {
            sss = ms.at(sssName);
        }
    }

    /*!
     * Performs the implementation specific updates before the physics calculations.
     *
     * @param tst The timestep start and duration .
     */
    virtual void updateBefore(const TimestepTime& tst) = 0;
    /*!
     *  Performs the implementation specific updates after the physics calculations.
     *
     * @param tst the timestep start and duration. Note that in some sense this
     *            update occurs near the end of the timestep at time tst.start + tst.duration
     */
    virtual void updateAfter(const TimestepTime& tst) = 0;

    /*!
     * Merges the ice-ocean fluxes and ocean-atmosphere fluxes into a single field to be passed to a
     * slab-ocean implementation or an ocean model through a coupler.
     */
    void mergeFluxes(const TimestepTime& tst)
    {
        assert(m_couplingArrays.checkAllRegistered());

        auto& fwFlux = fwFluxAccessor.getAutoRW();
        auto& qNoSun = qNoSunAccessor.getAutoRW();
        auto& tauX = tauXAccessor.getAutoRW();
        auto& tauY = tauYAccessor.getAutoRW();
        auto& sFlux = sFluxAccessor.getAutoRW();
        auto& qswNet = qswNetAccessor.getAutoRW();
        const auto& tauXIO = tauXIOAccessor.getAutoRO();
        const auto& cice = ciceAccessor.getAutoRO();
        const auto& deltaHice = deltaHiceAccessor.getAutoRO();
        const auto& rain = rainAccessor.getAutoRO();
        const auto& sss = sssAccessor.getAutoRO();
        const auto& qswow = qswowAccessor.getAutoRO();
        const auto& tauYIO = tauYIOAccessor.getAutoRO();
        const auto& evap = evapAccessor.getAutoRO();
        const auto& qio = qioAccessor.getAutoRO();
        const auto& newIce = newIceAccessor.getAutoRO();
        const auto& tauYOW = tauYOWAccessor.getAutoRO();
        const auto& qow = qowAccessor.getAutoRO();
        const auto& qswBase = qswBaseAccessor.getAutoRO();
        const auto& deltaSmelt = deltaSmeltAccessor.getAutoRO();
        const auto& tauXOW = tauXOWAccessor.getAutoRO();

        const FloatType dt = tst.step.seconds();

        overElementsAuto(OVER_ELEMENTS_LAMBDA(const ElementIndex i) {
            // Heat fluxes - partitioned in solar and non-solar
            qswNet[i] = cice[i] * qswBase[i] + (1 - cice[i]) * qswow[i];
            qNoSun[i] = cice[i] * qio[i] + (1 - cice[i]) * qow[i] - qswNet[i];

            // Mass fluxes - fresh water and salt
            // ice volume change, both laterally and vertically
            const FloatType deltaIceVol = newIce[i] + deltaHice[i];
            // the device compiler does not like a global constant appearing in the argument list of
            // a template function: "identifier "Ice::s" is undefined in device code"
            const FloatType s = Ice::s;
            // Effective ice salinity is always less than or equal to the SSS, and here we use
            // the right units too
            const FloatType effectiveIceSal = 1e-3 * Utils::min(s, sss[i]);

            // Positive flux is up!
            fwFlux[i]
                = ((1 - effectiveIceSal) * Ice::rho * deltaIceVol + Ice::rhoSnow * deltaSmelt[i])
                    / dt
                + (evap[i] - rain[i]) * (1 - cice[i]);
            sFlux[i] = effectiveIceSal * Ice::rho * deltaIceVol / dt;

            // Momentum fluxes
            tauX[i] = cice[i] * tauXIO[i] + (1 - cice[i]) * tauXOW[i];
            tauY[i] = cice[i] * tauYIO[i] + (1 - cice[i]) * tauYOW[i];
        });
    }

private:
protected:
    ModelArrayStore m_couplingArrays;

    ModelArrayAccessor<Shared::Q_IO, RW> qioAccessor; // Ice-ocean heat flux, W m⁻²
    ModelArrayAccessor<Protected::SST, RW>
        sstAccessor; // Coupled or slab ocean sea surface temperature, ˚C
    ModelArrayAccessor<Protected::SSS, RW>
        sssAccessor; // Coupled or slab ocean sea surface salinity, PSU
    ModelArrayAccessor<Protected::MLD, RW> mldAccessor; // Mixed layer or slab ocean depth m
    ModelArrayAccessor<Protected::TF, RW> tfAccessor; // Freezing point of the mixed layer, ˚C
    ModelArrayAccessor<Protected::ML_BULK_CP, RW>
        cpmlAccessor; // Heat capacity of the mixed layer, J K⁻¹ m²
    ModelArrayAccessor<Protected::OCEAN_U, RW> uAccessor; // x(east)-ward ocean current, m s⁻¹
    ModelArrayAccessor<Protected::OCEAN_V, RW> vAccessor; // y(north)-ward ocean current, m s⁻¹
    ModelArrayAccessor<Protected::SSH, RW> sshAccessor; // sea surface height, m
    ModelArrayAccessor<CouplingFields::Q_SS_NO_SW, RW>
        qNoSunAccessor; // Net surface ocean heat flux, except short wave, W m⁻²
    ModelArrayAccessor<CouplingFields::Q_SS_SW, RW>
        qswNetAccessor; // Net surface ocean shortwave flux, W m⁻²
    ModelArrayAccessor<CouplingFields::FWFLUX, RW>
        fwFluxAccessor; // Net surface ocean fresh-water flux, kg m⁻²
    ModelArrayAccessor<CouplingFields::SFLUX, RW>
        sFluxAccessor; // Net surface ocean salt flux, kg m⁻²
    ModelArrayAccessor<Shared::Q_SW_OW, RW> qswowAccessor; // Shortwave flux in open water W m⁻²
    ModelArrayAccessor<Shared::Q_SW_BASE, RW>
        qswBaseAccessor; // Shortwave flux at the base of the ice W m⁻²
    ModelArrayAccessor<CouplingFields::O_STRESS_X, RW>
        tauXAccessor; // x(east)-ward total ocean stress, Pa
    ModelArrayAccessor<CouplingFields::O_STRESS_Y, RW>
        tauYAccessor; // y(north)-ward total ocean stress, Pa

    ModelArrayAccessor<Shared::C_ICE_DG, RO> ciceAccessor;
    ModelArrayAccessor<Protected::IO_STRESS_X> tauXIOAccessor;
    ModelArrayAccessor<Protected::IO_STRESS_X> tauYIOAccessor;
    ModelArrayAccessor<Shared::EVAP, RW> evapAccessor;
    ModelArrayAccessor<Shared::RAIN, RO> rainAccessor;
    ModelArrayAccessor<Shared::NEW_ICE, RW> newIceAccessor;
    ModelArrayAccessor<Shared::DELTA_HICE, RW> deltaHiceAccessor;
    ModelArrayAccessor<Shared::HSNOW_MELT, RW> deltaSmeltAccessor;
    ModelArrayAccessor<Shared::Q_OW, RW> qowAccessor;
    ModelArrayAccessor<Shared::OW_STRESS_X, RW> tauXOWAccessor;
    ModelArrayAccessor<Shared::OW_STRESS_Y, RW> tauYOWAccessor;
};
} /* namespace Nextsim */

#endif /* IOCEANBOUNDARY_HPP */
