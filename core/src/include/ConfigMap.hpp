/*!
 * @file ConfigMap.hpp
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
typedef std::map<std::string, Fusi> ConfigMap;

enum ConfigMapType {
    DOUBLE,
    UNSIGNED,
    INT,
    STRING,
    N_TYPES
};
}
#endif /* CONFIGMAP_HPP */
