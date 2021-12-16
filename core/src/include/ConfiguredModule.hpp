/*!
 * @file ConfiguredModule.hpp
 *
 * @date Oct 29, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_CONFIGUREDMODULE_HPP
#define SRC_INCLUDE_CONFIGUREDMODULE_HPP

#include <string>

namespace Nextsim {

/*!
 * @brief A helper class to automatically load modules based on the
 * configuration.
 */
class ConfiguredModule {
public:
    ConfiguredModule();
    virtual ~ConfiguredModule();

    //! Parse the configuration for all of the modules defined in ModuleLoader.
    static void parseConfigurator();

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

#endif /* SRC_INCLUDE_CONFIGUREDMODULE_HPP */
