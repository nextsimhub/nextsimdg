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

    void update(const TimestepTime& tst);

private:
    ModelArrayRef<ProtectedArray::T_ICE> tice;
    ModelArrayRef<ProtectedArray::HTRUE_SNOW> h_snow_true; // cell-averaged value
};
}

#endif // NEXTSIM_DG_MONTHLYFLUXES_HPP
