/*!
 * @file CommandLineParser.cpp
 *
 * @date Oct 8, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/CommandLineParser.hpp"

#include <boost/program_options.hpp>
#include <cstdlib>
#include <iostream>
#include <string>

namespace Nextsim {

const std::string CommandLineParser::allModuleString = "all";
const std::string CommandLineParser::availableModuleString = "avail";

/*!
 * Parse the command line for command line options.
 *
 * @param argc the count of passed arguments.
 * @param argv The array of C strings passed as arguments.
 */
CommandLineParser::CommandLineParser(int argc, char* argv[])
{
    char oneConfigFile[] = "config-file";
    char manyConfigFiles[] = "config-files";
    char helpConfigStr[] = "help-config";

    boost::program_options::options_description opt("neXtSIM_DG command line options:");
    // clang-format off
    std::string helpConfigHelp = std::string("arg = [") + availableModuleString + " | moduleName | " + allModuleString + "]\n"
            "print help associated with one or more modules.\n"
            "avail    list all available module names.\n"
            "moduleName    print the help message associated with the given module.\n"
            "all    print all help messages associated with all available modules.";
    opt.add_options()
            ("help,h", "print help message")
            (helpConfigStr, boost::program_options::value<std::string>(), helpConfigHelp.c_str())
            (oneConfigFile, boost::program_options::value<std::string>(),
                    "specify a configuration file")
            (manyConfigFiles,
                    boost::program_options::value<std::vector<std::string>>()->multitoken(),
                    "specify a list of configuration files" )
             ;
    // clang-format on
    auto parsed = boost::program_options::command_line_parser(argc, argv)
                      .options(opt)
                      .style(boost::program_options::command_line_style::
                              unix_style) // | po::command_line_style::allow_long_disguise)
                      .allow_unregistered()
                      .run();
    boost::program_options::store(parsed, m_arguments);

    // Print help and exit
    if (m_arguments.count("help") || m_arguments.empty()) {
        std::cerr << opt << std::endl;
        std::exit(EXIT_SUCCESS);
    }

    m_configHelp.clear();

    // Store the config file names in order. The variables_map does not.
    for (auto& lionOption : parsed.options) {
        if (lionOption.string_key == oneConfigFile || lionOption.string_key == manyConfigFiles) {
            for (auto& stringValue : lionOption.value) {
                m_configFilenames.push_back(stringValue);
            }
        } else if (lionOption.string_key == helpConfigStr) {
            m_configHelp = lionOption.value.front();
        }
    }
}

/*!
 * Return a std::vector of the file names declared as config files on the command line
 */
std::vector<std::string> CommandLineParser::getConfigFileNames() { return m_configFilenames; }

} /* namespace Nextsim */
