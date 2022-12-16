/*!
 * @file Configurator.cpp
 *
 * @brief The source file for the Configurator class.
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
static NoAdditionalConfiguration noAddConf;
Configurator::AdditionalConfiguration* Configurator::p_addConf = &noAddConf;

// Parses the command line and named configuration file streams and returns the
// boost::program_options::variables_map.
boost::program_options::variables_map Configurator::parse(
    const boost::program_options::options_description& opt)
{
    boost::program_options::variables_map vm;

    // Parse the command file for any overrides
    int use_argc;
    char** use_argv;
    if (m_argc && m_argv) {
        use_argc = m_argc;
        use_argv = m_argv;
    } else {
        // Use a fake command line to ensure at least one parse happens in all cases
        use_argc = 1;
        char name[] = { 'n', 'e', 'x', 't', 's', 'i', 'm', 'd', 'g', '\0' };
        char* fake_argv[] = { name, nullptr };
        use_argv = fake_argv;
    }
    auto parsed = boost::program_options::command_line_parser(use_argc, use_argv)
                      .options(opt)
                      .style(boost::program_options::command_line_style::
                              unix_style) // | po::command_line_style::allow_long_disguise)
                      .allow_unregistered()
                      .run();
    boost::program_options::store(parsed, vm);

    // Parse the named streams for configuration
    for (auto iter = sources.begin(); iter != sources.end(); ++iter) {
        try {
            boost::program_options::store(
                boost::program_options::parse_config_file(**iter, opt, true), vm);
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

// Adds a string stream to the configuration sources. Useful for configuring
// classes in tests.
void Configurator::addSStream(const std::stringstream& sstream)
{
    addStream(std::move(std::unique_ptr<std::istream>(new std::stringstream(sstream.str()))));
}

// Sets the object providing additional configuration.
void Configurator::setAdditionalConfiguration(AdditionalConfiguration* pAC) { p_addConf = pAC; }

// Gets additional configuration from any object that exists.
void Configurator::getAdditionalConfiguration(const std::string& source)
{
    if (p_addConf) {
        addSStream(p_addConf->read(source));
    }
}

} /* namespace Nextsim */
