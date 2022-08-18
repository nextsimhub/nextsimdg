/*!
 * @file main.cpp
 * @date 11 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include <iostream>

#include "include/CommandLineParser.hpp"
#include "include/ConfigurationHelpPrinter.hpp"
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

    // Parse the configuration to load those that are explicitly configured
    Nextsim::ConfiguredModule::parseConfigurator();

    if (!cmdLine.configHelp().empty()) {
        Nextsim::Model::HelpMap map;
        Nextsim::Model::getHelpRecursive(map, true);
        Nextsim::ConfigurationHelpPrinter::setOutput(Nextsim::ConfigurationHelpPrinter::Output::ANSI);
        Nextsim::ConfigurationHelpPrinter::print(std::cout, map);
    } else {
        // Construct the Model
        Nextsim::Model model;
        // Apply the model configuration
        model.configure();
        // Run the Model
        model.run();
    }

    return 0;
}
