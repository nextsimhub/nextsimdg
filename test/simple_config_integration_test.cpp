/*!
 * @file simple_cfg_integration_test.cpp
 *
 * @date Nov 15, 2022
 * @author Alexander Smith <as3402@cam.ac.uk>
 */

#include "ArgV.hpp"
#include "CommandLineParser.hpp"
#include "ConfiguredModule.hpp"
#include "Configurator.hpp"
#include "ConfigMap.hpp"
#include "Model.hpp"

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <string>
#include <filesystem>

TEST_CASE("Read configuration file integration test", 
            "[Integration, CommandLineParser, Configurator]")
{

    const std::string logFilePrefix = "nextsim.";
    const std::string logFileSuffix = ".log";
 
    const std::string diagnosticPrefix = "diagnostic.";
    const std::string ncSuffix = ".nc";   

    const std::string restartFile = "restart.nc";
                
    const std::string config_files[2] 
        = {"config_simple_example.cfg", "config_simple_example_with_outdir.cfg"};

    const std::string initWorkingDir = std::filesystem::current_path();


    for (std::string targetConfigFilename: config_files) {
    
        DYNAMIC_SECTION("For config file " << targetConfigFilename << ": are parameters stored " 
                        << "in model and used correctly")
        {
            // Open example config file and load target answers from file 
            std::ifstream targetConfig(targetConfigFilename);
            std::string junk, value;

            std::vector<std::string> values;

            // skip header row
            targetConfig >> junk;

            // store values from last column only
            while (targetConfig >> junk >> junk >> value) {
                values.push_back(value);
            }

            const std::string initFileString = values[0];
            const std::string startString = values[1];
            const std::string stopString = values[2];
            const std::string timestepString = values[3];

            std::string outdirString;
            if (values.size() == 5) {
                outdirString = values[4];
                REQUIRE( outdirString.size() > 0 );
            } else {
                outdirString = "";
            }

            std::filesystem::path targetDir = std::filesystem::current_path();
            targetDir /= outdirString;

            // Require setup completed --- that the config_file for this test could be read correctly
            // before evaluating test sections
            REQUIRE( initFileString.size() > 0 );
            REQUIRE( startString.size() > 0 );
            REQUIRE( stopString.size() > 0 );
            REQUIRE( timestepString.size() > 0 );

            // Setup the command line arguments using this test class
            Nextsim::ArgV argv({"nextsimdg", "--config-file", targetConfigFilename});

            // Setup the configurator
            Nextsim::Configurator::setCommandLine(argv.argc(), argv());
            Nextsim::CommandLineParser cmdLine(argv.argc(), argv());
            Nextsim::Configurator::addFiles(cmdLine.getConfigFileNames());
            Nextsim::ConfiguredModule::parseConfigurator();

            // Setup the nextsim model
            Nextsim::Model model;
            model.configure();
            
            // Test that the config map stored in the model matches expected values as stored in the
            // config file
            
            SECTION( "Validate mapping of config map stored on model" )
            {
                Nextsim::ConfigMap cfgMap = model.getConfig();

                // Test strings stored in model config map match input config file
                // For Model::restartOptionName -> 'init_file'      
                std::string modelCfgInitFileString =  std::get<3>(cfgMap.find(model.keyMap.find( 
                                                        model.RESTARTFILE_KEY)->second)->second);
                REQUIRE_THAT(modelCfgInitFileString, Catch::Matchers::Equals(initFileString));

                // For 'start'
                std::string modelCfgStartString =  std::get<3>(cfgMap.find(model.keyMap.find( 
                                                        model.STARTTIME_KEY)->second)->second);
                REQUIRE_THAT(modelCfgStartString, Catch::Matchers::Equals(startString));

                // For 'stop'
                std::string modelCfgStopString =  std::get<3>(cfgMap.find(model.keyMap.find( 
                                                        model.STOPTIME_KEY)->second)->second);
                REQUIRE_THAT(modelCfgStopString, Catch::Matchers::Equals(stopString));

                // For 'time_step'
                std::string modelCfgTimestepString =  std::get<3>(cfgMap.find(model.keyMap.find( 
                                                        model.TIMESTEP_KEY)->second)->second);
                REQUIRE_THAT(modelCfgTimestepString, Catch::Matchers::Equals(timestepString));

                // For 'output_dir'
                std::string modelCfgOutdirString =  std::get<3>(cfgMap.find(model.keyMap.find( 
                                                        model.OUTPUTDIR_KEY)->second)->second);
                REQUIRE_THAT(modelCfgOutdirString, Catch::Matchers::Equals(outdirString));
            }


            SECTION("Run model with simple config")
            {
                model.run();

                // Check the target output dir exists (if specified in cfg, ensure dir was created)
                REQUIRE(std::filesystem::is_directory(targetDir));

                // Check output files were created in the right location for diagnostic, restart, log
                // Check the restart file exists
                REQUIRE(std::filesystem::exists(targetDir / restartFile));

                // Check the diagnostic file was generated for the last (only) step 
                std::string diagnosticFile = diagnosticPrefix + stopString + ncSuffix;
                REQUIRE(std::filesystem::exists(targetDir / diagnosticFile));

                // Check the nextsim.*.log file was generated. Instead of trying to determine the 
                // timestamp from when the model was run, loop over contents to ensure files
                // with the right prefix and extension were created.
                int logCounter = 0;
                for(auto &filename:std::filesystem::directory_iterator(targetDir)) {
                   if (filename.path().string().find(logFilePrefix) != std::string::npos &&
                      filename.path().string().find(logFileSuffix) != std::string::npos) {
                        logCounter += 1;
                        break;
                   }
                }
                REQUIRE(logCounter > 0);
            }
        

            //Teardown for dynamic section

            //If write location has been created/exists
            if (std::filesystem::is_directory(targetDir)) {

                //If the location, when resolved, is within the initial working directory
                if (std::filesystem::canonical(targetDir).string().find(initWorkingDir) == 0) {
                
                    for(auto &filename:std::filesystem::directory_iterator(targetDir)) {
                        //Remove log file    
                        if (filename.path().string().find(logFilePrefix) != std::string::npos &&
                              filename.path().string().find(logFileSuffix) != std::string::npos) {
                              std::filesystem::remove(filename.path());
                        }
                        
                        //Remove diagnostic file
                        if (filename.path().string().find(diagnosticPrefix) != std::string::npos &&
                              filename.path().string().find(ncSuffix) != std::string::npos) {
                               std::filesystem::remove(filename.path());  
                        }

                        //Remove restart file
                        if (filename.path().string().find(restartFile) != std::string::npos && 
                               filename.path().string().find(ncSuffix) != std::string::npos)  {
                               std::filesystem::remove(filename.path());
                        }
                    }
                }
            }
        } //end dynamic section
    }
}//end test case    
