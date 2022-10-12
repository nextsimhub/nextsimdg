//
// Created by Einar Ólason on 01/09/2022.
//

#ifndef NEXTSIM_DG_MU71ATMOSPHERE_HPP
#define NEXTSIM_DG_MU71ATMOSPHERE_HPP

#include "include/IIceAlbedo.hpp"
#include "include/constants.hpp"
#include "include/IAtmosphereBoundary.hpp"

namespace Nextsim {

class MU71Atmosphere : public IAtmosphereBoundary, public Configured<MU71Atmosphere> {

public:
    MU71Atmosphere()
        : tice(getProtectedArray()),
        h_snow_true(getProtectedArray())
    {};

    /*!
     * The required update call for an IFluxCalculation implementation. Here we just call
     * calculateElement inside an overElements loop.
     * @param tst The TimestepTime object for the current time step
     */
    void update(const TimestepTime& tst) override;

    void configure() override;

private:
    /*!
     * A function to calculate the fluxes qia, qio, and qow (=0), and subl (=0), as well as dqia_dt.
     * All incoming fluxes are tabulated from Maykut and Untersteiner (1971). Outgoing long wave and
     * the derivative (dqio_dt) are black body radiation.
     * @param i The index of the current grid cell
     * @param tst The TimestepTime object for the current time step
     */
    void calculateElement(size_t i, const TimestepTime& tst);

    ModelArrayRef<ProtectedArray::T_ICE, MARConstBackingStore> tice;
    ModelArrayRef<ProtectedArray::HTRUE_SNOW, MARConstBackingStore>
        h_snow_true; // cell-averaged value

    IIceAlbedo* iIceAlbedoImpl;
    double snowfall(const double dayOfYear, bool isLeap);
};

}

#endif // NEXTSIM_DG_MU71ATMOSPHERE_HPP
