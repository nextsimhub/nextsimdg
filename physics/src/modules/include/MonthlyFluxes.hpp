//
// Created by Einar Ólason on 01/09/2022.
//

#ifndef NEXTSIM_DG_MONTHLYFLUXES_HPP
#define NEXTSIM_DG_MONTHLYFLUXES_HPP

#include "include/IFluxCalculation.hpp"
#include "include/IIceAlbedo.hpp"
#include "include/constants.hpp"

namespace Nextsim {

class MonthlyFluxes : public IFluxCalculation {

public:
    MonthlyFluxes()
        : IFluxCalculation()
    {
    }

    /*!
     * The required update call for an IFluxCalculation implementation. Here we just call
     * calculateElement inside an overElements loop.
     * @param tst The TimestepTime object for the current time step
     */
    void update(const TimestepTime& tst);

private:
    /*!
     * A function to calculate the fluxes qia, qio, and qow (=0), and subl (=0), as well as dqia_dt.
     * All incoming fluxes are tabulated from Maykut and Untersteiner (1971). Outgoing long wave and
     * the derivative (dqio_dt) are black body radiation.
     * @param i The index of the current grid cell
     * @param tst The TimestepTime object for the current time step
     */
    void calculateElement(size_t i, const TimestepTime& tst);

    // We need acces to tice and h_snow_true
    ModelArrayRef<ProtectedArray::T_ICE> tice;
    ModelArrayRef<ProtectedArray::HTRUE_SNOW> h_snow_true; // cell-averaged value
};

}

#endif // NEXTSIM_DG_MONTHLYFLUXES_HPP
