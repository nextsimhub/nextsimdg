/*!
 * @file CommandLineParser.hpp
 *
 * @date Oct 8, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_COMMANDLINEPARSER_HPP
#define SRC_INCLUDE_COMMANDLINEPARSER_HPP

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

private:
    CommandLineParser() = default;

    boost::program_options::variables_map m_arguments;
    std::vector<std::string> m_configFilenames;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_COMMANDLINEPARSER_HPP */
