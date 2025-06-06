/*!
 * @file ILateralIceSpread.hpp
 *
 * @date 04 Jun 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#ifndef ILATERALICESPREAD_HPP
#define ILATERALICESPREAD_HPP

#include "include/ModelArrayRef.hpp"
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
        deltaCIce.resize();
        newice.resize();
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
        : deltaCIce(ModelArray::Type::H)
        , newice(ModelArray::Type::H)
        , snowMelt(ModelArray::Type::H)
        , cice(getStore())
        , deltaHi(getStore())
        , hice(getStore())
        , hsnow(getStore())
        , qow(getStore())
    {
        getStore().registerArray(Shared::DELTA_CICE, &deltaCIce, RW);
        getStore().registerArray(Shared::HSNOW_MELT, &snowMelt, RW);
        getStore().registerArray(Shared::NEW_ICE, &newice, RW);
    }

    HField deltaCIce; // Change in ice concentration
    HField newice; // New ice over open water this timestep, m
    HField snowMelt; // Ocean to snow transfer of freshwater kg m⁻²

    ModelArrayRef<Shared::C_ICE, RW> cice; // From IceGrowth
    ModelArrayRef<Shared::DELTA_HICE, RO> deltaHi; // From Vertical Ice Growth
    ModelArrayRef<Shared::H_ICE, RW> hice; // From IceGrowth
    ModelArrayRef<Shared::H_SNOW, RW> hsnow; // From Ice Growth?
    ModelArrayRef<Shared::Q_OW, RW> qow; // From FluxCalculation
};

} /* namespace Nextsim */

#endif /* ILATERALICESPREAD_HPP */
