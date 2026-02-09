/*
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/FiniteElementSpecHum.hpp"

#include <cmath>

namespace Nextsim {

FiniteElementSpecHum FiniteElementSpecHum::m_water(
    6.1121e2, 18.729, 257.87, 227.3, 7.2e-4, 3.20e-6, 5.9e-10);
FiniteElementSpecHum FiniteElementSpecHum::m_ice(
    6.1115e2, 23.036, 279.82, 333.7, 2.2e-4, 3.83e-6, 6.4e-10);

} /* namespace Nextsim */
