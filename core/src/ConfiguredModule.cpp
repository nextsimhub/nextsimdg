/*!
 * @file ConfiguredModule.cpp
 *
 * @date Oct 29, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ConfiguredModule.hpp"

#include "include/Configurator.hpp"
#include "include/ModuleLoader.hpp"

#include <boost/program_options.hpp>
#include <stdexcept>
namespace Nextsim {

const std::string ConfiguredModule::MODULE_PREFIX = "Modules";

ConfiguredModule::ConfiguredModule()
{
    // TODO Auto-generated constructor stub
}

ConfiguredModule::~ConfiguredModule()
{
    // TODO Auto-generated destructor stub
}

void ConfiguredModule::parseConfigurator()
{
    // A default string that can never be a valid C++ class name
    std::string defaultStr = "+++DEFAULT+++";
    // Construct a new options map
    boost::program_options::options_description opt;

    ModuleLoader& loader = ModuleLoader::getLoader();
    for (const std::string& module : loader.listModules()) {
        std::string defaultImpl = *loader.listImplementations(module).begin();
        opt.add_options()(addPrefix(module).c_str(),
            boost::program_options::value<std::string>()->default_value(defaultStr),
            ("Load an implementation of " + module).c_str());
    }

    boost::program_options::variables_map vm = Configurator::parse(opt);

    for (const std::string& module : loader.listModules()) {
        std::string implString = vm[addPrefix(module)].as<std::string>();
        // Only do anything if the retrieved option is not the default value
        if (implString != defaultStr) {
            // Check that the retrieved value is one of the implementations
            // defined for this module
            std::string moduleImpl = "";
            for (const std::string implName : loader.listImplementations(module)) {
                if (implString == implName) moduleImpl = implName;
            }
            if (moduleImpl != "") {
                loader.setImplementation(module, implString);
            } else {
                std::string what = "Invalid implementation \"";
                what += implString + "\" of module " + module +".";
                throw std::domain_error(what);
            }
        }
    }
}

std::string ConfiguredModule::addPrefix(const std::string& moduleName)
{
    return MODULE_PREFIX + "." + moduleName;
}

} /* namespace Nextsim */
