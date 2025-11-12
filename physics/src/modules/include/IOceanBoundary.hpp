/*!
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IOCEANBOUNDARY_HPP
#define IOCEANBOUNDARY_HPP

#include "include/CheckingModelComponent.hpp"
#include "include/ModelArrayAccessor.hpp"
#include "include/constants.hpp"
#include "include/gridNames.hpp"

namespace Nextsim {

//! An interface class for the oceanic inputs into the ice physics.
class IOceanBoundary : public CheckingModelComponent {
public:
    IOceanBoundary()
        : qioAccessor(getStore(), RW, ModelArray::Type::H, std::pair(-1e8, 1e8))
        , sstAccessor(getStore(), RO, ModelArray::Type::H, std::pair(-5.0, 50.0))
        , sssAccessor(getStore(), RO, ModelArray::Type::H, std::pair(0.0, 50.0))
        , mldAccessor(getStore(), RO, ModelArray::Type::H, std::pair(1e-3, 12e3))
        , cpmlAccessor(getStore(), RO, ModelArray::Type::H, std::pair(0.0, 1e11))
        , tfAccessor(getStore(), RO, ModelArray::Type::H, std::pair(-5.0, 0.0))
        , uAccessor(getStore(), RO, ModelArray::Type::H, std::pair(-1.0, 1.0))
        , vAccessor(getStore(), RO, ModelArray::Type::H, std::pair(-1.0, 1.0))
        , sshAccessor(getStore(), RO, ModelArray::Type::H, std::pair(-10.0, 10.0))
        , qNoSunAccessor(m_couplingArrays, RO, ModelArray::Type::H, std::pair(-1e6, 1e6))
        , qswNetAccessor(m_couplingArrays, RO, ModelArray::Type::H, std::pair(-1e3, 1e3))
        , fwFluxAccessor(m_couplingArrays, RO, ModelArray::Type::H)
        , sFluxAccessor(m_couplingArrays, RO, ModelArray::Type::H)
        , qswowAccessor(getStore(), RW, ModelArray::Type::H, std::pair(-1e3, 1e-6))
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
        qio.resize();
        HField& sst = sstAccessor.getHostRW();
        sst.resize();
        HField& sss = sssAccessor.getHostRW();
        sss.resize();
        HField& mld = mldAccessor.getHostRW();
        mld.resize();
        HField& cpml = cpmlAccessor.getHostRW();
        cpml.resize();
        HField& tf = tfAccessor.getHostRW();
        tf.resize();
        UField& u = uAccessor.getHostRW();
        u.resize();
        VField& v = vAccessor.getHostRW();
        v.resize();
        HField& ssh = sshAccessor.getHostRW();
        ssh.resize();
        HField& qNoSun = qNoSunAccessor.getHostRW();
        qNoSun.resize();
        HField& qswNet = qswNetAccessor.getHostRW();
        qswNet.resize();
        HField& fwFlux = fwFluxAccessor.getHostRW();
        fwFlux.resize();
        HField& sFlux = sFluxAccessor.getHostRW();
        sFlux.resize();
        HField& qswow = qswowAccessor.getHostRW();
        qswow.resize();
        HField& qswBase = qswBaseAccessor.getHostRW();
        qswBase.resize();
        HField& tauX = tauXAccessor.getHostRW();
        tauX.resize();
        HField& tauY = tauYAccessor.getHostRW();
        tauY.resize();

        addChecks({
            { "qio", &qio },
            { "sst", &sst },
            { "sss", &sss },
            { "mld", &mld },
            { "cpml", &cpml },
            { "tf", &tf },
            { "u", &u },
            { "v", &v },
            { "ssh", &ssh },
            { "qNoSun", &qNoSun },
            { "qswNet", &qswNet },
            { "fwFlux", &fwFlux },
            { "sFlux", &sFlux },
            { "qswow", &qswow },
            { "qswBase", &qswBase },
            { "tauX", &tauX },
            { "tauY", &tauY },
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

        const double dt = tst.step.seconds();
        HField& fwFlux = fwFluxAccessor.getHostRW();
        HField& qNoSun = qNoSunAccessor.getHostRW();
        HField& tauX = tauXAccessor.getHostRW();
        HField& tauY = tauYAccessor.getHostRW();
        HField& sFlux = sFluxAccessor.getHostRW();
        HField& qswNet = qswNetAccessor.getHostRW();
        const HField& tauXIO = tauXIOAccessor.getHostRO();
        const HField& cice = ciceAccessor.getHostRO();
        const HField& deltaHice = deltaHiceAccessor.getHostRO();
        const HField& rain = rainAccessor.getHostRO();
        const HField& sss = sssAccessor.getHostRO();
        const HField& qswow = qswowAccessor.getHostRO();
        const HField& tauYIO = tauYIOAccessor.getHostRO();
        const HField& evap = evapAccessor.getHostRO();
        const HField& qio = qioAccessor.getHostRO();
        const HField& newIce = newIceAccessor.getHostRO();
        const HField& tauYOW = tauYOWAccessor.getHostRO();
        const HField& qow = qowAccessor.getHostRO();
        const HField& qswBase = qswBaseAccessor.getHostRO();
        const HField& deltaSmelt = deltaSmeltAccessor.getHostRO();
        const HField& tauXOW = tauXOWAccessor.getHostRO();

        overElements(
            [&](const size_t i, const TimestepTime& tsTime) {
                // Heat fluxes - partitioned in solar and non-solar
                qswNet[i] = cice[i] * qswBase[i] + (1 - cice[i]) * qswow[i];
                qNoSun[i] = cice[i] * qio[i] + (1 - cice[i]) * qow[i] - qswNet[i];

                // Mass fluxes - fresh water and salt
                // ice volume change, both laterally and vertically
                const double deltaIceVol = newIce[i] + deltaHice[i] * cice[i];
                // change in snow volume due to melting (should be < 0)
                const double meltSnowVol = deltaSmelt[i] * cice[i];
                // Effective ice salinity is always less than or equal to the SSS, and here we use
                // the right units too
                const double effectiveIceSal = 1e-3 * std::min(Ice::s, sss[i]);

                // Positive flux is up!
                fwFlux[i]
                    = ((1 - effectiveIceSal) * Ice::rho * deltaIceVol + Ice::rhoSnow * meltSnowVol)
                        / dt
                    + (evap[i] - rain[i]) * (1 - cice[i]);
                sFlux[i] = effectiveIceSal * Ice::rho * deltaIceVol / dt;

                // Momentum fluxes
                tauX[i] = cice[i] * tauXIO[i] + (1 - cice[i]) * tauXOW[i];
                tauY[i] = cice[i] * tauYIO[i] + (1 - cice[i]) * tauYOW[i];
            },
            tst);
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
