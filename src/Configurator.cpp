/*!
 * @file Configurator.cpp
 *
 * @date Oct 8, 2021
 * @author Tim Spain
 */

#include "include/Configurator.hpp"

#include <iostream>

namespace Nextsim {

std::vector<std::unique_ptr<std::istream>> Configurator::sources;
int Configurator::m_argc;
char** Configurator::m_argv;

boost::program_options::variables_map Configurator::parse(const boost::program_options::options_description& opt)
{
    boost::program_options::variables_map vm;

    // Parse the command file for any overrides
    if (m_argc && m_argv) {
        auto parsed = boost::program_options::command_line_parser(m_argc, m_argv)
                                          .options(opt)
                                          .style(boost::program_options::command_line_style::unix_style)// | po::command_line_style::allow_long_disguise)
                                          .allow_unregistered()
                                          .run();
        boost::program_options::store(parsed, vm);
    }

    // Parse the named streams for configuration
    for (auto iter = sources.begin(); iter != sources.end(); ++iter) {
        try {
        boost::program_options::store( boost::program_options::parse_config_file(**iter, opt, true), vm);
        } catch (std::exception& e) {
            // Echo the exception, but carry on
            std::cerr << e.what() << std::endl;
        }
        // Once the stream has been parsed, clear all flags (especially EOF)
        // and seek back to the start.
        (*iter)->clear();
        (*iter)->seekg(0, std::ios_base::beg);
    }

    return vm;
}
} /* namespace Nextsim */
