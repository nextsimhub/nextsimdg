/*!
 * @file ConfiguredModule.hpp
 *
 * @date Oct 29, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CONFIGUREDMODULE_HPP
#define CONFIGUREDMODULE_HPP

#include "include/ConfigMap.hpp"

#include <functional>
#include <list>
#include <map>
#include <string>
#include <utility>

namespace Nextsim {

/*!
 * @brief A helper class to load modules based on the configuration.
 */
class ConfiguredModule {
public:
    ConfiguredModule() = default;
    virtual ~ConfiguredModule() = default;

    /*!
     * Gets the implementation name for a named module.
     * Returns an empty string if there is no such module.
     *
     * @param interface the name of the module to be checked.
     */
    static std::string getImpl(const std::string& interface);

    /*!
     * Returns a map from interface name to implementation name for all modules
     * configured in this class.
     *
     */
    static ConfigMap getAllModuleConfigurations();
    /*!
     * @brief Transforms a module name into a configuration option name.
     *
     * @param moduleName The module name to be prefixed.
     */
    static std::string addPrefix(const std::string& moduleName);

    //! The configuration options section name for modules.
    static const std::string MODULE_PREFIX;
};

} /* namespace Nextsim */

#endif /* CONFIGUREDMODULE_HPP */
