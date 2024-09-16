/*!
 * @file ConfiguredModule.cpp
 *
 * @date Oct 29, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ConfiguredModule.hpp"

#include "include/Configurator.hpp"
#include "include/NextsimModule.hpp"

#include <boost/program_options.hpp>
#include <stdexcept>

namespace Nextsim {

const std::string ConfiguredModule::MODULE_PREFIX = "Modules";

// A default string that can never be a valid C++ class name
static const std::string defaultStr = "+++DEFAULT+++";

std::string ConfiguredModule::addPrefix(const std::string& moduleName)
{
    return MODULE_PREFIX + "." + moduleName;
}

std::string ConfiguredModule::getImpl(const std::string& module)
{
    boost::program_options::options_description opt;

    opt.add_options()(addPrefix(module).c_str(),
        boost::program_options::value<std::string>()->default_value(defaultStr),
        ("Load an implementation of " + module).c_str());

    std::string implString = Configurator::parse(opt)[addPrefix(module)].as<std::string>();
    // Only do anything if the retrieved option is not the default value
    if (implString != defaultStr) {
        return implString;
    } else {
        return "";
    }
}

ConfigMap ConfiguredModule::getAllModuleConfigurations()
{
    ConfigMap iiMap;
    for (const auto& moduleImpl : Module::ImplementationNames::getAll()) {
        iiMap[moduleImpl.first] = moduleImpl.second;
    }
    return iiMap;
}
} /* namespace Nextsim */
