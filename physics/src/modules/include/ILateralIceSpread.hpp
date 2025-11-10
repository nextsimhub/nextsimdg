/*!
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Einar Ólason <einar.olason@nersc.no>
 */

#ifndef ILATERALICESPREAD_HPP
#define ILATERALICESPREAD_HPP

#include "include/ModelArrayAccessor.hpp"
#include "include/ModelComponent.hpp"
#include "include/Time.hpp"

namespace Nextsim {
//! An interface class to calculate the freezing of open water and melting of ice.
class ILateralIceSpread : public ModelComponent {
public:
    virtual ~ILateralIceSpread() = default;

    std::string getName() const override { return "LateralIceSpread"; }
    void setData(const ModelState::DataMap& ms) override
    {
        HField& deltaCIce = deltaCIceAccessor.getHostRW();
        deltaCIce.resize();
        HField& newice = newiceAccessor.getHostRW();
        newice.resize();
        HField& snowMelt = snowMeltAccessor.getHostRW();
        snowMelt.resize();

        /*
         * Set these to zero, so we don't have uninitialized values floating around.
         * As long as we have the nextsim_thermo.use_thermo_forcing option, all
         * derived classes must initialise these values to something sensible.
         */
        deltaCIce = 0.;
        newice = 0.;
        snowMelt = 0.;
    }

    /*!
     * Updates the lateral ice melt and formation for the timestep.
     *
     * @param tStep The object containing the timestep start and duration times.
     */
    virtual void update(const TimestepTime& tstep) = 0;

protected:
    ILateralIceSpread()
        : deltaCIceAccessor(getStore(), RW, ModelArray::Type::H)
        , newiceAccessor(getStore(), RW, ModelArray::Type::H)
        , snowMeltAccessor(getStore(), RW, ModelArray::Type::H)
        , ciceAccessor(getStore())
        , deltaHiAccessor(getStore())
        , hiceAccessor(getStore())
        , hsnowAccessor(getStore())
        , qowAccessor(getStore())
    {
    }

    ModelArrayAccessor<Shared::DELTA_CICE, RW> deltaCIceAccessor; // Change in ice concentration
    ModelArrayAccessor<Shared::NEW_ICE, RW>
        newiceAccessor; // New ice over open water this timestep, m
    ModelArrayAccessor<Shared::HSNOW_MELT, RW>
        snowMeltAccessor; // Ocean to snow transfer of freshwater kg m⁻²

    ModelArrayAccessor<Shared::C_ICE_DG, RW> ciceAccessor; // From IceGrowth
    ModelArrayAccessor<Shared::DELTA_HICE, RO> deltaHiAccessor; // From Vertical Ice Growth
    ModelArrayAccessor<Shared::H_ICE_DG, RW> hiceAccessor; // From IceGrowth
    ModelArrayAccessor<Shared::H_SNOW_DG, RW> hsnowAccessor; // From Ice Growth?
    ModelArrayAccessor<Shared::Q_OW, RW> qowAccessor; // From FluxCalculation
};

} /* namespace Nextsim */

#endif /* ILATERALICESPREAD_HPP */
