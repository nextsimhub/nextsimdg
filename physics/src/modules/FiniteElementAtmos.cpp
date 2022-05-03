/*
 * FiniteElementAtmos.cpp
 *
 *  Created on: 2 May 2022
 *      Author: work
 */

#include "include/FiniteElementAtmos.hpp"

#include "include/constants.hpp"
#include "include/ISpecificHumidity.hpp"

namespace Nextsim {
void FiniteElementAtmos::update(const TimestepTime& tst)
{
    overElements(
        std::bind(&FiniteElementAtmos::updateElement, this, std::placeholders::_1, std::placeholders::_2),
        tst);
}

void FiniteElementAtmos::updateElement(size_t i, const TimestepTime& tst)
{
    // Specific heat of...
    // ...the air
    sh_air[i] = SpecificHumidity::water()(t_dew2[i], p_air[i]);
    // ...over the open ocean
    sh_water[i] = SpecificHumidity::water()(sst[i], p_air[i], sss[i]);
    // ...over the ice
    sh_ice[i] = SpecificHumidity::ice()(tice.zIndexAndLayer(i, 0), p_air[i]);

    // Density of the wet air
    double Ra_wet = Air::Ra / (1 - sh_air[i] * (1 - Vapour::Ra / Air::Ra));
    rho_air[i] = p_air[i] / (Ra_wet * kelvin(t_air[i]));

    // Heat capacity of the wet air
    cp_wet[i] = Air::cp + sh_air[i] * Vapour::cp;
}

} /* namespace Nextsim */
