/*!
 * @file ConfiguredModule.cpp
 *
 * @date Oct 29, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ConfiguredModule.hpp"

#include "include/Configurator.hpp"
#include "include/Module.hpp"

#include <boost/program_options.hpp>
#include <stdexcept>

namespace Nextsim {

const std::string ConfiguredModule::MODULE_PREFIX = "Modules";

ConfiguredModule::map ConfiguredModule::configuredModules;

// A default string that can never be a valid C++ class name
static const std::string defaultStr = "+++DEFAULT+++";

void ConfiguredModule::parseConfigurator()
{
    // Construct a new options map
    boost::program_options::options_description opt;

    for (auto entry : configuredModules) {
        std::string module = entry.first;
        opt.add_options()(addPrefix(module).c_str(),
            boost::program_options::value<std::string>()->default_value(defaultStr),
            ("Load an implementation of " + module).c_str());
    }

    boost::program_options::variables_map vm = Configurator::parse(opt);

    for (auto entry : configuredModules) {
        std::string implString = vm[addPrefix(entry.first)].as<std::string>();
        // Only do anything if the retrieved option is not the default value
        if (implString != defaultStr) {
            entry.second.first(implString);
        }
    }
}

std::string ConfiguredModule::addPrefix(const std::string& moduleName)
{
    return MODULE_PREFIX + "." + moduleName;
}

void ConfiguredModule::setConfiguredModules(const map& ls)
{
    map newMap(ls);
    configuredModules.merge(newMap);
}

void ConfiguredModule::configureModule(const std::string& mod, const fn& setter, const ofn& getter)
{
    configuredModules[mod] = { setter, getter };
}

std::string ConfiguredModule::getImpl(const std::string& module)
{
    if (cachedImplementations().count(module))
        return cachedImplementations().at(module);
    else {
        boost::program_options::options_description opt;

        opt.add_options()(addPrefix(module).c_str(),
                boost::program_options::value<std::string>()->default_value(defaultStr),
                ("Load an implementation of " + module).c_str());

        std::string implString = Configurator::parse(opt)[addPrefix(module)].as<std::string>();
        // Only do anything if the retrieved option is not the default value
        if (implString != defaultStr) {
            cachedImplementations()[module] = implString;
            return implString;
        } else {
            return "";
        }
    }
}

ConfigMap ConfiguredModule::getAllModuleConfigurations()
{
    ConfigMap iiMap;
    for (auto entry : configuredModules) {
        iiMap[addPrefix(entry.first)] = entry.second.second();
    }
    return iiMap;
}
} /* namespace Nextsim */
