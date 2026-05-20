/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/IceMinima.hpp"

namespace Nextsim {

const FloatType IceMinima::hMinDefault = 0.01;
const FloatType IceMinima::cMinDefault = 1e-12;

FloatType IceMinima::hMin = hMinDefault;
FloatType IceMinima::cMin = cMinDefault;

} /* namespace Nextsim */
