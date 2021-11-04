/*!
 * @file main.cpp
 * @date 11 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include <iostream>

#include "include/CommandLineParser.hpp"
#include "include/Configurator.hpp"
#include "include/ConfiguredModule.hpp"
#include "include/Model.hpp"

int main(int argc, char* argv[])
{

    // Pass the command line to Configurator to handle
    Nextsim::Configurator::setCommandLine(argc, argv);
    // Extract any config files defined on the command line
    Nextsim::CommandLineParser cmdLine(argc, argv);
    // Pass the config file names to Configurator
    Nextsim::Configurator::addFiles(cmdLine.getConfigFileNames());

    // Translate the configuration to the ModuleLoader
    Nextsim::ConfiguredModule::parseConfigurator();

    // Construct the Model
    Nextsim::Model model = Nextsim::Model();
    // Run the Model
    model.run();

    return 0;
}
