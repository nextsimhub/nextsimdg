/*!
 * @file IOceanBoundary.hpp
 *
 * @date 29 Apr 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IOCEANBOUNDARY_HPP
#define IOCEANBOUNDARY_HPP

#include "include/ModelComponent.hpp"
#include "include/constants.hpp"

namespace Nextsim {

//! An interface class for the oceanic inputs into the ice physics.
class IOceanBoundary : public ModelComponent {
public:
    IOceanBoundary()
        : qio(ModelArray::Type::H)
        , sst(ModelArray::Type::H)
        , sss(ModelArray::Type::H)
        , mld(ModelArray::Type::H)
        , cpml(ModelArray::Type::H)
        , tf(ModelArray::Type::H)
        , u(ModelArray::Type::H)
        , v(ModelArray::Type::H)
        , ssh(ModelArray::Type::H)
        , qNoSun(ModelArray::Type::H)
        , qswNet(ModelArray::Type::H)
        , fwFlux(ModelArray::Type::H)
        , sFlux(ModelArray::Type::H)
        , qswow(ModelArray::Type::H)
        , qswBase(ModelArray::Type::H)
        , tauX(ModelArray::Type::H)
        , tauY(ModelArray::Type::H)
        , cice(getStore())
        , emp(getStore())
        , newIce(getStore())
        , deltaHice(getStore())
        , deltaSmelt(getStore())
        , qow(getStore())
        , tauXIO(getStore())
        , tauYIO(getStore())
        , tauXOW(getStore())
        , tauYOW(getStore())
    {
        // Receive
        m_couplingArrays.registerArray(CouplingFields::MLD, &mld, RW);
        m_couplingArrays.registerArray(CouplingFields::OCEAN_U, &u, RW);
        m_couplingArrays.registerArray(CouplingFields::OCEAN_V, &v, RW);
        m_couplingArrays.registerArray(CouplingFields::SSH, &ssh, RW);
        m_couplingArrays.registerArray(CouplingFields::SSS, &sss, RW);
        m_couplingArrays.registerArray(CouplingFields::SST, &sst, RW);
        // Send
        m_couplingArrays.registerArray(CouplingFields::FWFLUX, &fwFlux, RO);
        m_couplingArrays.registerArray(CouplingFields::O_STRESS_X, &tauX, RO);
        m_couplingArrays.registerArray(CouplingFields::O_STRESS_Y, &tauY, RO);
        m_couplingArrays.registerArray(CouplingFields::Q_SS_NO_SW, &qNoSun, RO);
        m_couplingArrays.registerArray(CouplingFields::Q_SS_SW, &qswNet, RO);
        m_couplingArrays.registerArray(CouplingFields::SFLUX, &sFlux, RO);

        getStore().registerArray(Shared::Q_IO, &qio, RW);
        getStore().registerArray(Protected::SST, &sst, RO);
        getStore().registerArray(Protected::SSS, &sss, RO);
        getStore().registerArray(Protected::MLD, &mld, RO);
        getStore().registerArray(Protected::ML_BULK_CP, &cpml, RO);
        getStore().registerArray(Protected::TF, &tf, RO);
        getStore().registerArray(Protected::OCEAN_U, &u, RO);
        getStore().registerArray(Protected::OCEAN_V, &v, RO);
        getStore().registerArray(Protected::SSH, &ssh, RO);
        getStore().registerArray(Shared::Q_SW_OW, &qswow, RW);
        getStore().registerArray(Shared::Q_SW_BASE, &qswBase, RW);
    }
    virtual ~IOceanBoundary() = default;

    ModelState getState() const override { return ModelState(); }
    ModelState getState(const OutputLevel&) const override { return getState(); }

    std::string getName() const override { return "IOceanBoundary"; }
    void setData(const ModelState::DataMap& ms) override
    {
        qio.resize();
        sst.resize();
        sss.resize();
        mld.resize();
        cpml.resize();
        tf.resize();
        u.resize();
        v.resize();
        ssh.resize();
        qNoSun.resize();
        qswNet.resize();
        fwFlux.resize();
        sFlux.resize();
        qswow.resize();
        qswBase.resize();
        tauX.resize();
        tauY.resize();

        if (ms.count("sst")) {
            sst = ms.at("sst");
        }
        if (ms.count("sss")) {
            sss = ms.at("sss");
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
        dt = tst.step.seconds();
        overElements(std::bind(&IOceanBoundary::mergeFluxesElement, this, std::placeholders::_1,
                         std::placeholders::_2),
            tst);
    }

private:
    double dt;

    void mergeFluxesElement(size_t i, const TimestepTime& tst)
    {
        // Heat fluxes - partitioned in solar and non-solar
        qswNet[i] = cice[i] * qswBase[i] + (1 - cice[i]) * qswow[i];
        qNoSun[i] = cice[i] * qio[i] + (1 - cice[i]) * qow[i] - qswNet[i];

        // Mass fluxes - fresh water and salt
        // ice volume change, both laterally and vertically
        const double deltaIceVol = newIce[i] + deltaHice[i] * cice[i];
        // change in snow volume due to melting (should be < 0)
        const double meltSnowVol = deltaSmelt[i] * cice[i];
        // Effective ice salinity is always less than or equal to the SSS, and here we use the right
        // units too
        const double effectiveIceSal = 1e-3 * std::min(Ice::s, sss[i]);

        // Positive flux is up!
        fwFlux[i]
            = ((1 - effectiveIceSal) * Ice::rho * deltaIceVol + Ice::rhoSnow * meltSnowVol) / dt
            + emp[i] * (1 - cice[i]);
        sFlux[i] = effectiveIceSal * Ice::rho * deltaIceVol / dt;

        // Momentum fluxes
        tauX[i] = cice[i] * tauXIO[i] + (1 - cice[i]) * tauXOW[i];
        tauY[i] = cice[i] * tauYIO[i] + (1 - cice[i]) * tauYOW[i];
    }

protected:
    HField qio; // Ice-ocean heat flux, W m⁻²
    HField sst; // Coupled or slab ocean sea surface temperature, ˚C
    HField sss; // Coupled or slab ocean sea surface salinity, PSU
    HField mld; // Mixed layer or slab ocean depth m
    HField tf; // Freezing point of the mixed layer, ˚C
    HField cpml; // Heat capacity of the mixed layer, J K⁻¹ m²
    UField u; // x(east)-ward ocean current, m s⁻¹
    VField v; // y(north)-ward ocean current, m s⁻¹
    HField ssh; // sea surface height, m
    HField qNoSun; // Net surface ocean heat flux, except short wave, W m⁻²
    HField qswNet; // Net surface ocean shortwave flux, W m⁻²
    HField fwFlux; // Net surface ocean fresh-water flux, kg m⁻²
    HField sFlux; // Net surface ocean salt flux, kg m⁻²
    HField qswow; // Shortwave flux in open water W m⁻²
    HField qswBase; // Shortwave flux at the base of the ice W m⁻²
    HField tauX; // x(east)-ward total ocean stress, Pa
    HField tauY; // y(north)-ward total ocean stress, Pa

    ModelArrayReferenceStore m_couplingArrays;

    ModelArrayRef<Protected::C_ICE, RO> cice;
    ModelArrayRef<Protected::EVAP_MINUS_PRECIP, RO> emp;
    ModelArrayRef<Protected::IO_STRESS_X> tauXIO;
    ModelArrayRef<Protected::IO_STRESS_X> tauYIO;
    ModelArrayRef<Shared::NEW_ICE, RW> newIce;
    ModelArrayRef<Shared::DELTA_HICE, RW> deltaHice;
    ModelArrayRef<Shared::HSNOW_MELT, RW> deltaSmelt;
    ModelArrayRef<Shared::Q_OW, RW> qow;
    ModelArrayRef<Shared::OW_STRESS_X, RW> tauXOW;
    ModelArrayRef<Shared::OW_STRESS_Y, RW> tauYOW;
};
} /* namespace Nextsim */

#endif /* IOCEANBOUNDARY_HPP */
