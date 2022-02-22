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

void ConfiguredModule::parseConfigurator()
{
    // A default string that can never be a valid C++ class name
    std::string defaultStr = "+++DEFAULT+++";
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
            entry.second(implString);
        }
    }
}

std::string ConfiguredModule::addPrefix(const std::string& moduleName)
{
    return MODULE_PREFIX + "." + moduleName;
}

void ConfiguredModule::setConfiguredModules(const map& ls)
{
    for (auto entry : ls) {
        configureModule(entry.first, entry.second);
    }
}

void ConfiguredModule::configureModule(const std::string& mod, const fn& function)
{
    configuredModules[mod] = function;
}
} /* namespace Nextsim */
