/*!
 * @file gridNames.hpp
 *
 * @date Oct 24, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef GRIDNAMES_HPP
#define GRIDNAMES_HPP

#include <string>

namespace Nextsim {

// A set of canonical names for fields used in file reading and writing.
static const std::string hiceName = "hice";
static const std::string ciceName = "cice";
static const std::string hsnowName = "hsnow";
static const std::string ticeName = "tice";
static const std::string maskName = "mask";

static const std::string coordsName = "coords";

static const std::string mdiName = "missing_value";

static const std::string timeName = "time";

}

#endif /* GRIDNAMES_HPP */
