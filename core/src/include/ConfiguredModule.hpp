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
 * @brief A helper class to automatically load modules based on the
 * configuration.
 */
class ConfiguredModule {
public:
    typedef std::function<void(const std::string&)> fn;
    typedef std::function<std::string()> ofn;
    typedef std::map<const std::string, std::pair<fn, ofn>> map;
    ConfiguredModule() = default;
    virtual ~ConfiguredModule() = default;

    //! Parse the configuration for all of the modules defined in ModuleLoader.
    static void parseConfigurator();

    /*!
     * Sets the names of all the modules to be configured
     * @param ls a map of Module names and the locations of implementation
     *           setter and getter functions.
     */
    static void setConfiguredModules(const map& ls);

    /*!
     * Adds the name of a module to the list to be configured
     * @param mod the name of the class to be configured.
     * @param setter the void(std::string) function to set the implementation.
     * @param getter the std::string() function to get the implementation.
     */
    static void configureModule(const std::string& mod, const fn& setter, const ofn& getter);

    /*!
     * Gets the implementation name for a named module.
     * Returns an empty string if there is no such module.
     *
     * @param interface the name of the module to be checked.
     */
    static std::string getModuleConfiguration(const std::string& interface);

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

private:
    static map configuredModules;
};

} /* namespace Nextsim */

#endif /* CONFIGUREDMODULE_HPP */
