/*!
 * @file main.cpp
 * @date 11 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include <iostream>
#include <mpi.h>
#include <oasis_c.h>


#include "include/CommandLineParser.hpp"
#include "include/ConfigurationHelpPrinter.hpp"
#include "include/Configurator.hpp"
#include "include/ConfiguredModule.hpp"
#include "include/Model.hpp"
#include "include/NetcdfMetadataConfiguration.hpp"

int main(int argc, char* argv[])
{
    // OASIS init
    char *comp_name = "nextsim";
    int comp_id;
    OASIS_CHECK_ERR(oasis_c_init_comp(&comp_id, comp_name, OASIS_COUPLED));

    // OASIS definition
    // OASIS defining partition
    int part_params[OASIS_Serial_Params];
    int part_in;
    part_params[OASIS_Strategy] = OASIS_Serial;
    part_params[OASIS_Offset] = 0;
    part_params[OASIS_Length] = 30*30;
    OASIS_CHECK_ERR(oasis_c_def_partition(&part_in, OASIS_Serial_Params,
        				  part_params, part_params[OASIS_Length],
					  "part_in"));
    // OASIS defining variable
    int bundle_size;
    int variable;
    bundle_size = 1;
    OASIS_CHECK_ERR(oasis_c_def_var(&variable, "FRECVNXT", part_in, bundle_size,
  				    OASIS_IN, OASIS_REAL));
    // OASIS finalising definition
    OASIS_CHECK_ERR(oasis_c_enddef());

    // Pass the command line to Configurator to handle
    Nextsim::Configurator::setCommandLine(argc, argv);
    // Extract any config files defined on the command line
    Nextsim::CommandLineParser cmdLine(argc, argv);
    // Pass the config file names to Configurator
    Nextsim::Configurator::addFiles(cmdLine.getConfigFileNames());
    try {
        // Get the configuration stored in the restart file
        Nextsim::NetcdfMetadataConfiguration ncdfMC;
        Nextsim::Configurator::setAdditionalConfiguration(&ncdfMC);
        Nextsim::Configurator::getAdditionalConfiguration(Nextsim::Model::restartOptionName);
    } catch (const std::exception& e) {
        // Do nothing. If there is no additional configuration to be parse, ignore it.
    }
    // Parse the configuration to load those that are explicitly configured
    Nextsim::ConfiguredModule::parseConfigurator();

    if (!cmdLine.configHelp().empty()) {
        Nextsim::Model::HelpMap map;
        Nextsim::Model::getHelpRecursive(map, true);
        Nextsim::ConfigurationHelpPrinter::setOutput(
            Nextsim::ConfigurationHelpPrinter::Output::ANSI);
        Nextsim::ConfigurationHelpPrinter::print(std::cout, map, cmdLine.configHelp());
    } else {
        // Construct the Model
        Nextsim::Model model;
        // Apply the model configuration
        model.configure();
        // Run the Model
        model.run();
    }
    // OASIS receiving field 
    int kinfo;
    int date=0;
    double bundle[30][30][1];
    OASIS_CHECK_ERR(oasis_c_get(variable, date, 30, 30, 1, OASIS_REAL, OASIS_ROW_MAJOR, bundle, &kinfo));

    // OASIS ending simulation
    OASIS_CHECK_ERR(oasis_c_terminate());
    return 0;
}
