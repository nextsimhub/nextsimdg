/*!
 * @file ModelState.hpp
 *
 * @date Feb 28, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_INCLUDE_MODELSTATE_HPP
#define CORE_SRC_INCLUDE_MODELSTATE_HPP

#include <map>
#include <string>
#include "ModelArray.hpp"

namespace Nextsim {

typedef std::map<std::string, const ModelArray*> ModelState;

} /* namespace Nextsim */

#endif /* CORE_SRC_INCLUDE_MODELSTATE_HPP */
