/*!
 * @file ConfiguredModule.hpp
 *
 * @date Oct 29, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_CONFIGUREDMODULE_HPP
#define SRC_INCLUDE_CONFIGUREDMODULE_HPP

#include <functional>
#include <list>
#include <map>
#include <string>

namespace Nextsim {

/*!
 * @brief A helper class to automatically load modules based on the
 * configuration.
 */
class ConfiguredModule {
public:

    typedef std::function<void(const std::string&)> fn;
    typedef std::map<const std::string, fn> map;
    ConfiguredModule() = default;
    virtual ~ConfiguredModule() = default;

    //! Parse the configuration for all of the modules defined in ModuleLoader.
    static void parseConfigurator();

    /*!
     * Sets the names of all the modules to be configured
     * @param nameList a std::list of std::strings with the names of the
     *                  interfaces to be configured.
     */
    static void setConfiguredModules(const map& ls);

    /*!
     * Adds the name of a module to the list to be configured
     * @param name the name of the class to be configured.
     */
    static void configureModule(const std::string& mod, const fn&);

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

#endif /* SRC_INCLUDE_CONFIGUREDMODULE_HPP */
