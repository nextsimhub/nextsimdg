/*!
 * @file ConfigMap.hpp
 *
 * @brief A header file defining the ConfigMap type.
 *
 * @date Aug 19, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CONFIGMAP_HPP
#define CONFIGMAP_HPP

#include <map>
#include <string>
#include <variant>

namespace Nextsim {

typedef std::variant<double, unsigned, int, std::string> Fusi;
/*!
 * A typedef type that allows model configuration items to be stored.
 */
typedef std::map<std::string, Fusi> ConfigMap;

static const std::size_t CONFIGMAP_DOUBLE = 0;
static const std::size_t CONFIGMAP_UNSIGNED = 1;
static const std::size_t CONFIGMAP_INT = 2;
static const std::size_t CONFIGMAP_STRING = 3;

}
#endif /* CONFIGMAP_HPP */
