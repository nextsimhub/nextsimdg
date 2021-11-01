/*
 * @file ConfiguredModule.cpp
 *
 * @date Oct 29, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ConfiguredModule.hpp"

#include "include/Configurator.hpp"
#include "include/ModuleLoader.hpp"

#include <boost/program_options.hpp>
namespace Nextsim {

ConfiguredModule::ConfiguredModule() {
  // TODO Auto-generated constructor stub
}

ConfiguredModule::~ConfiguredModule() {
  // TODO Auto-generated destructor stub
}

void ConfiguredModule::parseConfigurator()
{
    // Construct a new options map
    boost::program_options::options_description opt;

    ModuleLoader& loader = ModuleLoader::getLoader();
    for (const std::string& module : loader.listModules()) {
        std::string defaultImpl = *loader.listImplementations(module).begin();
        opt.add_options()
        (module.c_str(), boost::program_options::value<std::string>()->default_value(defaultImpl), ("Load an implementation of " + module).c_str())
        ;
    }

    boost::program_options::variables_map vm = Configurator::parse(opt);

    for (const std::string& module : loader.listModules()) {
        std::string implString = vm[module].as<std::string>();
        loader.setImplementation(module, implString);
    }
}
} /* namespace Nextsim */
