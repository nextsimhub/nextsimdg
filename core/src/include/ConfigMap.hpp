/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CONFIGMAP_HPP
#define CONFIGMAP_HPP

#include <map>
#include <string>
#include <variant>

namespace Nextsim {

typedef std::variant<double, float, unsigned, int, std::string> Fusi;
typedef std::map<std::string, Fusi> ConfigMap;

// clang-format off
enum ConfigMapType {
    DOUBLE,
    FLOAT,
    UNSIGNED,
    INT,
    STRING,
    N_TYPES
};
// clang-format on
}
#endif /* CONFIGMAP_HPP */
