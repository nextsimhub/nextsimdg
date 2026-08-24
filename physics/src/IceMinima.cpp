/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/IceMinima.hpp"

namespace Nextsim {

const FloatType IceMinima::hMinDefault = 0.01;
const FloatType IceMinima::cMinDefault = std::is_same_v<FloatType, float> ? 1e-5 : 1e-12;

FloatType IceMinima::hMin = hMinDefault;
FloatType IceMinima::cMin = cMinDefault;

} /* namespace Nextsim */
