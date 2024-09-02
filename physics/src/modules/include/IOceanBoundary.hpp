/*!
 * @file IOceanBoundary.hpp
 *
 * @date 30 Aug 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IOCEANBOUNDARY_HPP
#define IOCEANBOUNDARY_HPP

#include "include/ModelComponent.hpp"
#include "include/constants.hpp"

namespace Nextsim {

namespace CouplingFields {
    constexpr TextTag SST = "SST"; // sea surface temperature ˚C
    constexpr TextTag SSS = "SSS"; // sea surface salinity PSU
    constexpr TextTag MLD = "MLD"; // Mixed layer or slab ocean depth m
    constexpr TextTag OCEAN_U = "U"; // x(east)-ward ocean current m s⁻¹
    constexpr TextTag OCEAN_V = "V"; // y(north)-ward ocean current m s⁻¹
}
//! An interface class for the oceanic inputs into the ice physics.
class IOceanBoundary : public ModelComponent {
public:
    IOceanBoundary()
        : cice(getStore())
        , qswow(getStore())
        , emp(getStore())
        , newIce(getStore())
        , deltaHice(getStore())
        , deltaSmelt(getStore())
        , qow(getStore())
        , qswBase(getStore())
    {
        m_couplingArrays.registerArray(CouplingFields::SST, &sst, RW);
        m_couplingArrays.registerArray(CouplingFields::SSS, &sss, RW);
        m_couplingArrays.registerArray(CouplingFields::OCEAN_U, &u, RW);
        m_couplingArrays.registerArray(CouplingFields::OCEAN_V, &v, RW);

        getStore().registerArray(Shared::Q_IO, &qio, RW);
        getStore().registerArray(Protected::SST, &sst, RO);
        getStore().registerArray(Protected::SSS, &sss, RO);
        getStore().registerArray(Protected::MLD, &mld, RO);
        getStore().registerArray(Protected::ML_BULK_CP, &cpml, RO);
        getStore().registerArray(Protected::TF, &tf, RO);
        getStore().registerArray(Protected::OCEAN_U, &u, RO);
        getStore().registerArray(Protected::OCEAN_V, &v, RO);
        getStore().registerArray(Protected::FWFLUX, &fwFlux, RO);
        getStore().registerArray(Protected::SFLUX, &sFlux, RO);
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
        qNoSun.resize();
        qswNet.resize();
        fwFlux.resize();
        sFlux.resize();

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
    void mergeFluxes(size_t i, const TimestepTime& tst)
    {
        const double dt = tst.step.seconds();

        qswNet(i) = cice(i) * qswBase(i) + (1 - cice(i)) * qswow(i);
        qNoSun(i) = cice(i) * qio(i) + (1 - cice(i)) * qow(i) - qswNet(i);

        // ice volume change, both laterally and vertically
        const double deltaIceVol = newIce(i) + deltaHice(i) * cice(i);
        // change in snow volume due to melting (should be < 0)
        const double meltSnowVol = deltaSmelt(i) * cice(i);
        // Effective ice salinity is always less than or equal to the SSS, and here we use the right
        // units too
        const double effectiveIceSal = 1e-3 * std::min(sss(i), Ice::s);

        // Positive flux is up!
        fwFlux(i)
            = ((1 - effectiveIceSal) * Ice::rho * deltaIceVol + Ice::rhoSnow * meltSnowVol) / dt
            + emp(i) * (1 - cice(i));
        sFlux(i) = effectiveIceSal * Ice::rho * deltaIceVol / dt;
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
    HField qNoSun; // Net surface ocean heat flux, except short wave, W m⁻²
    HField qswNet; // Net surface ocean shortwave flux, W m⁻²
    HField fwFlux; // Net surface ocean fresh-water flux, kg m⁻²
    HField sFlux; // Net surface ocean salt flux, kg m⁻²

    ModelArrayReferenceStore m_couplingArrays;

    ModelArrayRef<Protected::C_ICE, RO> cice;
    ModelArrayRef<Shared::Q_SW_OW, RO> qswow;
    ModelArrayRef<Protected::EVAP_MINUS_PRECIP, RO> emp;
    ModelArrayRef<Shared::NEW_ICE, RW> newIce;
    ModelArrayRef<Shared::DELTA_HICE, RW> deltaHice;
    ModelArrayRef<Shared::HSNOW_MELT, RW> deltaSmelt;
    ModelArrayRef<Shared::Q_OW, RW> qow;
    ModelArrayRef<Shared::Q_SW_BASE, RW> qswBase;
};
} /* namespace Nextsim */

#endif /* IOCEANBOUNDARY_HPP */
