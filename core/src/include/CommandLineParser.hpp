/*!
 * @file CommandLineParser.hpp
 *
 * @date Oct 8, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef COMMANDLINEPARSER_HPP
#define COMMANDLINEPARSER_HPP

#include <boost/program_options.hpp>
#include <string>
#include <vector>

namespace Nextsim {

/*!
 * A class to parse the neXtSIM_DG command line.
 */
class CommandLineParser {
public:
    /*!
     * Parse the command line for command line options.
     *
     * @param argc the count of passed arguments.
     * @param argv The array of C strings passed as arguments.
     */
    CommandLineParser(int argc, char* argv[]);
    virtual ~CommandLineParser() = default;

    /*!
     * Return a std::vector of the file names declared as config files on the
     * command line, in order.
     */
    std::vector<std::string> getConfigFileNames();

    /*!
     * Returns the name of any module for which help was requested.
     */
    std::string configHelp() { return m_configHelp; }

    /*!
     * The special string denoting help for all available modules.
     */
    const static std::string allModuleString;

    /*!
     * The special string denoting a request for listing available modules.
     */
    const static std::string availableModuleString;

private:
    CommandLineParser() = default;

    boost::program_options::variables_map m_arguments;
    std::vector<std::string> m_configFilenames;
    std::string m_configHelp;
};

} /* namespace Nextsim */

#endif /* COMMANDLINEPARSER_HPP */
